//! Shared pipeline cache format (`kira-organelle.bin`) per `CACHE_FILE.md`.

use std::fs::File;
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use memmap2::Mmap;
use sprs::CsMat;
use thiserror::Error;
use tracing::info;

use crate::input::InputError;
use crate::io::mtx::{MtxInput, load_mtx_dir};

const HEADER_SIZE: usize = 256;
const ALIGNMENT: usize = 64;
const MAGIC: &[u8; 4] = b"KORG";
const VERSION_MAJOR: u16 = 1;
const VERSION_MINOR: u16 = 0;
const ENDIAN_TAG: u32 = 0x1234_5678;

/// Errors for shared organelle cache files.
#[derive(Debug, Error)]
pub enum OrganelleCacheError {
    #[error("I/O error reading {path:?}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("cache format error in {path:?}: {message}")]
    Format { path: PathBuf, message: String },
    #[error("input parsing error: {0}")]
    Input(#[from] InputError),
}

/// Read-only mmap-backed view of an organelle cache file.
#[derive(Debug)]
pub struct OrganelleCacheView {
    mmap: Mmap,
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    col_ptr_offset: usize,
    row_idx_offset: usize,
    values_u32_offset: usize,
}

impl OrganelleCacheView {
    /// CSC `col_ptr`.
    pub fn col_ptr(&self) -> &[u64] {
        let len = self.n_cells + 1;
        let bytes = &self.mmap[self.col_ptr_offset..self.col_ptr_offset + len * 8];
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u64, len) }
    }

    /// CSC `row_idx`.
    pub fn row_idx(&self) -> &[u32] {
        let bytes = &self.mmap[self.row_idx_offset..self.row_idx_offset + self.nnz * 4];
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u32, self.nnz) }
    }

    /// CSC `values_u32`.
    pub fn values_u32(&self) -> &[u32] {
        let bytes = &self.mmap[self.values_u32_offset..self.values_u32_offset + self.nnz * 4];
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u32, self.nnz) }
    }

    /// Convert the view to an in-memory MTX-like structure used by the pipeline.
    pub fn to_mtx_input(&self) -> MtxInput {
        let indptr: Vec<usize> = self.col_ptr().iter().map(|v| *v as usize).collect();
        let indices: Vec<usize> = self.row_idx().iter().map(|v| *v as usize).collect();
        let data: Vec<f32> = self.values_u32().iter().map(|v| *v as f32).collect();
        let matrix = CsMat::new_csc((self.n_genes, self.n_cells), indptr, indices, data);
        MtxInput {
            matrix,
            features: self.genes.clone(),
            barcodes: self.barcodes.clone(),
        }
    }
}

/// Build shared cache from MTX directory and write to destination path.
pub fn write_organelle_bin_from_mtx(
    out_path: &Path,
    input_path: &Path,
    gene_symbol_col: Option<usize>,
) -> Result<(), OrganelleCacheError> {
    let mtx = load_mtx_dir(input_path, gene_symbol_col)?;
    write_organelle_bin(out_path, &mtx)
}

/// Write shared cache from loaded MTX input.
pub fn write_organelle_bin(path: &Path, input: &MtxInput) -> Result<(), OrganelleCacheError> {
    info!(path = ?path, "Shared cache writer: start");
    let n_genes = input.matrix.rows();
    let n_cells = input.matrix.cols();
    let nnz = input.matrix.nnz();
    info!(
        n_genes,
        n_cells, nnz, "Shared cache writer: matrix dimensions"
    );

    if input.features.len() != n_genes {
        return Err(OrganelleCacheError::Format {
            path: path.to_path_buf(),
            message: "feature count does not match matrix rows".to_string(),
        });
    }
    if input.barcodes.len() != n_cells {
        return Err(OrganelleCacheError::Format {
            path: path.to_path_buf(),
            message: "barcode count does not match matrix cols".to_string(),
        });
    }

    // Enforce row ordering and value type for deterministic/portable CSC.
    for (col_idx, col) in input.matrix.outer_iterator().enumerate() {
        let mut prev_row: Option<usize> = None;
        for (row, value) in col.iter() {
            if let Some(prev) = prev_row
                && row <= prev
            {
                return Err(OrganelleCacheError::Format {
                    path: path.to_path_buf(),
                    message: format!("row_idx not strictly increasing in column {col_idx}"),
                });
            }
            prev_row = Some(row);
            if !value.is_finite()
                || *value < 0.0
                || value.fract() != 0.0
                || *value > u32::MAX as f32
            {
                return Err(OrganelleCacheError::Format {
                    path: path.to_path_buf(),
                    message: "values_u32 must be finite non-negative integer counts".to_string(),
                });
            }
        }
    }

    let genes_table = encode_string_table(&input.features)?;
    let barcodes_table = encode_string_table(&input.barcodes)?;
    info!(
        genes_table_bytes = genes_table.len(),
        barcodes_table_bytes = barcodes_table.len(),
        "Shared cache writer: encoded string tables"
    );

    let col_ptr_len = input.matrix.indptr().raw_storage().len();
    let row_idx_len = input.matrix.indices().len();
    let values_len = input.matrix.data().len();

    let genes_table_offset = align_to(HEADER_SIZE, ALIGNMENT);
    let barcodes_table_offset = align_to(genes_table_offset + genes_table.len(), ALIGNMENT);
    let col_ptr_offset = align_to(barcodes_table_offset + barcodes_table.len(), ALIGNMENT);
    let col_ptr_bytes = col_ptr_len * 8;
    let row_idx_offset = align_to(col_ptr_offset + col_ptr_bytes, ALIGNMENT);
    let row_idx_bytes = row_idx_len * 4;
    let values_offset = align_to(row_idx_offset + row_idx_bytes, ALIGNMENT);
    let values_bytes = values_len * 4;
    let file_bytes = values_offset + values_bytes;
    info!(
        file_bytes,
        col_ptr_offset, row_idx_offset, values_offset, "Shared cache writer: computed layout"
    );

    let mut header = [0u8; HEADER_SIZE];
    header[0..4].copy_from_slice(MAGIC);
    write_u16(&mut header, 4, VERSION_MAJOR);
    write_u16(&mut header, 6, VERSION_MINOR);
    write_u32(&mut header, 8, ENDIAN_TAG);
    write_u32(&mut header, 12, HEADER_SIZE as u32);
    write_u64(&mut header, 16, n_genes as u64);
    write_u64(&mut header, 24, n_cells as u64);
    write_u64(&mut header, 32, nnz as u64);
    write_u64(&mut header, 40, genes_table_offset as u64);
    write_u64(&mut header, 48, genes_table.len() as u64);
    write_u64(&mut header, 56, barcodes_table_offset as u64);
    write_u64(&mut header, 64, barcodes_table.len() as u64);
    write_u64(&mut header, 72, col_ptr_offset as u64);
    write_u64(&mut header, 80, row_idx_offset as u64);
    write_u64(&mut header, 88, values_offset as u64);
    write_u64(&mut header, 96, 0); // n_blocks
    write_u64(&mut header, 104, 0); // blocks_offset
    write_u64(&mut header, 112, file_bytes as u64);
    write_u64(&mut header, 128, 0); // data_crc64 (v1)
    let crc = crc64_ecma(&header);
    write_u64(&mut header, 120, crc);

    let mut file = File::create(path).map_err(|source| OrganelleCacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    file.set_len(file_bytes as u64)
        .map_err(|source| OrganelleCacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    write_at(&mut file, 0, &header, path)?;
    write_at(&mut file, genes_table_offset, &genes_table, path)?;
    write_at(&mut file, barcodes_table_offset, &barcodes_table, path)?;
    write_u64_from_usize_iter(
        &mut file,
        col_ptr_offset,
        input.matrix.indptr().raw_storage().iter().copied(),
        path,
    )?;
    write_u32_from_usize_iter(
        &mut file,
        row_idx_offset,
        input.matrix.indices().iter().copied(),
        path,
    )?;
    write_u32_from_f32_iter(
        &mut file,
        values_offset,
        input.matrix.data().iter().copied(),
        path,
    )?;
    file.sync_all().map_err(|source| OrganelleCacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    info!(path = ?path, "Shared cache writer: done");

    Ok(())
}

/// Memory-map and validate shared cache.
pub fn mmap_organelle_bin(path: &Path) -> Result<OrganelleCacheView, OrganelleCacheError> {
    let file = File::open(path).map_err(|source| OrganelleCacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let metadata = file.metadata().map_err(|source| OrganelleCacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() < HEADER_SIZE as u64 {
        return Err(format_error(path, "file too small"));
    }

    let mmap = unsafe {
        Mmap::map(&file).map_err(|source| OrganelleCacheError::Io {
            path: path.to_path_buf(),
            source,
        })?
    };
    let header = &mmap[..HEADER_SIZE];
    validate_header(path, header, metadata.len() as usize)?;

    let n_genes = read_u64(header, 16) as usize;
    let n_cells = read_u64(header, 24) as usize;
    let nnz = read_u64(header, 32) as usize;
    let genes_table_offset = read_u64(header, 40) as usize;
    let genes_table_bytes = read_u64(header, 48) as usize;
    let barcodes_table_offset = read_u64(header, 56) as usize;
    let barcodes_table_bytes = read_u64(header, 64) as usize;
    let col_ptr_offset = read_u64(header, 72) as usize;
    let row_idx_offset = read_u64(header, 80) as usize;
    let values_u32_offset = read_u64(header, 88) as usize;
    let file_bytes = read_u64(header, 112) as usize;

    bounds_check(path, genes_table_offset, genes_table_bytes, file_bytes)?;
    bounds_check(
        path,
        barcodes_table_offset,
        barcodes_table_bytes,
        file_bytes,
    )?;
    bounds_check(path, col_ptr_offset, (n_cells + 1) * 8, file_bytes)?;
    bounds_check(path, row_idx_offset, nnz * 4, file_bytes)?;
    bounds_check(path, values_u32_offset, nnz * 4, file_bytes)?;

    let genes = decode_string_table(
        &mmap[genes_table_offset..genes_table_offset + genes_table_bytes],
        n_genes,
        path,
    )?;
    let barcodes = decode_string_table(
        &mmap[barcodes_table_offset..barcodes_table_offset + barcodes_table_bytes],
        n_cells,
        path,
    )?;

    let view = OrganelleCacheView {
        mmap,
        n_genes,
        n_cells,
        nnz,
        genes,
        barcodes,
        col_ptr_offset,
        row_idx_offset,
        values_u32_offset,
    };
    validate_csc(path, &view)?;
    Ok(view)
}

fn validate_header(
    path: &Path,
    header: &[u8],
    actual_len: usize,
) -> Result<(), OrganelleCacheError> {
    if &header[0..4] != MAGIC {
        return Err(format_error(path, "invalid magic"));
    }
    if read_u16(header, 4) != VERSION_MAJOR {
        return Err(format_error(path, "unsupported major version"));
    }
    if read_u32(header, 8) != ENDIAN_TAG {
        return Err(format_error(path, "invalid endian tag"));
    }
    if read_u32(header, 12) as usize != HEADER_SIZE {
        return Err(format_error(path, "invalid header size"));
    }
    if read_u64(header, 96) != 0 || read_u64(header, 104) != 0 {
        return Err(format_error(
            path,
            "v1 requires n_blocks=0 and blocks_offset=0",
        ));
    }
    if read_u64(header, 128) != 0 {
        return Err(format_error(path, "v1 requires data_crc64=0"));
    }
    let file_bytes = read_u64(header, 112) as usize;
    if file_bytes != actual_len {
        return Err(format_error(path, "file_bytes does not match file size"));
    }
    let expected = read_u64(header, 120);
    let mut crc_bytes = [0u8; HEADER_SIZE];
    crc_bytes.copy_from_slice(header);
    write_u64(&mut crc_bytes, 120, 0);
    let actual = crc64_ecma(&crc_bytes);
    if expected != actual {
        return Err(format_error(path, "header_crc64 mismatch"));
    }
    Ok(())
}

fn validate_csc(path: &Path, view: &OrganelleCacheView) -> Result<(), OrganelleCacheError> {
    let col_ptr = view.col_ptr();
    let row_idx = view.row_idx();

    if col_ptr.is_empty() || col_ptr[0] != 0 {
        return Err(format_error(path, "col_ptr[0] must be 0"));
    }
    let nnz = view.nnz as u64;
    if *col_ptr.last().unwrap_or(&u64::MAX) != nnz {
        return Err(format_error(path, "col_ptr[n_cells] must equal nnz"));
    }
    for w in col_ptr.windows(2) {
        if w[0] > w[1] {
            return Err(format_error(
                path,
                "col_ptr must be monotonic non-decreasing",
            ));
        }
    }
    for col in 0..view.n_cells {
        let start = col_ptr[col] as usize;
        let end = col_ptr[col + 1] as usize;
        if end > row_idx.len() || start > end {
            return Err(format_error(path, "invalid col_ptr bounds"));
        }
        let mut prev: Option<u32> = None;
        for row in &row_idx[start..end] {
            if *row as usize >= view.n_genes {
                return Err(format_error(path, "row_idx out of bounds"));
            }
            if let Some(p) = prev
                && *row <= p
            {
                return Err(format_error(
                    path,
                    "row_idx must be strictly increasing inside each column",
                ));
            }
            prev = Some(*row);
        }
    }
    Ok(())
}

fn encode_string_table(values: &[String]) -> Result<Vec<u8>, OrganelleCacheError> {
    let count = values.len();
    let mut offsets = Vec::with_capacity(count + 1);
    offsets.push(0u32);
    let mut blob = Vec::new();
    for value in values {
        let bytes = value.as_bytes();
        if (blob.len() + bytes.len()) > u32::MAX as usize {
            return Err(OrganelleCacheError::Format {
                path: PathBuf::from("<memory>"),
                message: "string table too large".to_string(),
            });
        }
        blob.extend_from_slice(bytes);
        offsets.push(blob.len() as u32);
    }

    let mut out = Vec::with_capacity(4 + 4 * (count + 1) + blob.len());
    out.extend_from_slice(&(count as u32).to_le_bytes());
    for off in offsets {
        out.extend_from_slice(&off.to_le_bytes());
    }
    out.extend_from_slice(&blob);
    Ok(out)
}

fn decode_string_table(
    bytes: &[u8],
    expected_count: usize,
    path: &Path,
) -> Result<Vec<String>, OrganelleCacheError> {
    if bytes.len() < 4 {
        return Err(format_error(path, "string table too small"));
    }
    let count = read_u32(bytes, 0) as usize;
    if count != expected_count {
        return Err(format_error(path, "string table count mismatch"));
    }
    let offsets_bytes = 4 * (count + 1);
    if bytes.len() < 4 + offsets_bytes {
        return Err(format_error(path, "string table offsets out of bounds"));
    }
    let mut offsets = Vec::with_capacity(count + 1);
    for i in 0..=count {
        offsets.push(read_u32(bytes, 4 + i * 4) as usize);
    }
    for w in offsets.windows(2) {
        if w[0] > w[1] {
            return Err(format_error(path, "string offsets must be monotonic"));
        }
    }
    let blob = &bytes[4 + offsets_bytes..];
    if offsets[count] != blob.len() {
        return Err(format_error(path, "string offsets terminal value mismatch"));
    }
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let s = &blob[offsets[i]..offsets[i + 1]];
        let value = std::str::from_utf8(s).map_err(|_| format_error(path, "invalid UTF-8"))?;
        out.push(value.to_string());
    }
    Ok(out)
}

fn write_at(
    file: &mut File,
    offset: usize,
    bytes: &[u8],
    path: &Path,
) -> Result<(), OrganelleCacheError> {
    file.seek(SeekFrom::Start(offset as u64))
        .map_err(|source| OrganelleCacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    file.write_all(bytes)
        .map_err(|source| OrganelleCacheError::Io {
            path: path.to_path_buf(),
            source,
        })
}

fn bounds_check(
    path: &Path,
    offset: usize,
    len: usize,
    file_bytes: usize,
) -> Result<(), OrganelleCacheError> {
    let end = offset
        .checked_add(len)
        .ok_or_else(|| format_error(path, "section overflow"))?;
    if end > file_bytes {
        return Err(format_error(path, "section out of file bounds"));
    }
    Ok(())
}

fn format_error(path: &Path, message: &str) -> OrganelleCacheError {
    OrganelleCacheError::Format {
        path: path.to_path_buf(),
        message: message.to_string(),
    }
}

fn align_to(value: usize, alignment: usize) -> usize {
    let rem = value % alignment;
    if rem == 0 {
        value
    } else {
        value + (alignment - rem)
    }
}

fn write_u16(buf: &mut [u8], offset: usize, value: u16) {
    buf[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

fn write_u32(buf: &mut [u8], offset: usize, value: u32) {
    buf[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn write_u64(buf: &mut [u8], offset: usize, value: u64) {
    buf[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

fn read_u16(buf: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes(buf[offset..offset + 2].try_into().unwrap())
}

fn read_u32(buf: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(buf[offset..offset + 4].try_into().unwrap())
}

fn read_u64(buf: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap())
}

fn write_u64_from_usize_iter<I: Iterator<Item = usize>>(
    file: &mut File,
    offset: usize,
    iter: I,
    path: &Path,
) -> Result<(), OrganelleCacheError> {
    file.seek(SeekFrom::Start(offset as u64))
        .map_err(|source| OrganelleCacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    let mut buf = [0u8; 8 * 1024];
    let mut used = 0usize;
    for value in iter {
        let bytes = (value as u64).to_le_bytes();
        if used + 8 > buf.len() {
            file.write_all(&buf[..used])
                .map_err(|source| OrganelleCacheError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            used = 0;
        }
        buf[used..used + 8].copy_from_slice(&bytes);
        used += 8;
    }
    if used > 0 {
        file.write_all(&buf[..used])
            .map_err(|source| OrganelleCacheError::Io {
                path: path.to_path_buf(),
                source,
            })?;
    }
    Ok(())
}

fn write_u32_from_usize_iter<I: Iterator<Item = usize>>(
    file: &mut File,
    offset: usize,
    iter: I,
    path: &Path,
) -> Result<(), OrganelleCacheError> {
    file.seek(SeekFrom::Start(offset as u64))
        .map_err(|source| OrganelleCacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    let mut buf = [0u8; 4 * 1024];
    let mut used = 0usize;
    for value in iter {
        let bytes = (value as u32).to_le_bytes();
        if used + 4 > buf.len() {
            file.write_all(&buf[..used])
                .map_err(|source| OrganelleCacheError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            used = 0;
        }
        buf[used..used + 4].copy_from_slice(&bytes);
        used += 4;
    }
    if used > 0 {
        file.write_all(&buf[..used])
            .map_err(|source| OrganelleCacheError::Io {
                path: path.to_path_buf(),
                source,
            })?;
    }
    Ok(())
}

fn write_u32_from_f32_iter<I: Iterator<Item = f32>>(
    file: &mut File,
    offset: usize,
    iter: I,
    path: &Path,
) -> Result<(), OrganelleCacheError> {
    file.seek(SeekFrom::Start(offset as u64))
        .map_err(|source| OrganelleCacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    let mut buf = [0u8; 4 * 1024];
    let mut used = 0usize;
    for value in iter {
        let bytes = (value as u32).to_le_bytes();
        if used + 4 > buf.len() {
            file.write_all(&buf[..used])
                .map_err(|source| OrganelleCacheError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            used = 0;
        }
        buf[used..used + 4].copy_from_slice(&bytes);
        used += 4;
    }
    if used > 0 {
        file.write_all(&buf[..used])
            .map_err(|source| OrganelleCacheError::Io {
                path: path.to_path_buf(),
                source,
            })?;
    }
    Ok(())
}

fn crc64_ecma(bytes: &[u8]) -> u64 {
    let poly: u64 = 0x42F0_E1EB_A9EA_3693;
    let mut crc: u64 = 0;
    for byte in bytes {
        crc ^= (*byte as u64) << 56;
        for _ in 0..8 {
            if (crc & 0x8000_0000_0000_0000) != 0 {
                crc = (crc << 1) ^ poly;
            } else {
                crc <<= 1;
            }
        }
    }
    crc
}

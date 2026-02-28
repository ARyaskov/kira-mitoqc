//! Gene symbol column auto-detection for features.tsv.

use std::collections::HashSet;
use std::io::{BufRead, BufReader, Read};

use tracing::{error, info};

use crate::input::InputError;

struct ColumnStats {
    total: usize,
    mt_hits: usize,
    oxphos_hits: usize,
    ens_hits: usize,
    alpha_like_hits: usize,
    unique: HashSet<String>,
}

impl ColumnStats {
    fn new() -> Self {
        Self {
            total: 0,
            mt_hits: 0,
            oxphos_hits: 0,
            ens_hits: 0,
            alpha_like_hits: 0,
            unique: HashSet::new(),
        }
    }

    fn score(&self) -> f32 {
        if self.total == 0 {
            return 0.0;
        }
        let total = self.total as f32;
        let mt_rate = self.mt_hits as f32 / total;
        let ox_rate = self.oxphos_hits as f32 / total;
        let ens_rate = self.ens_hits as f32 / total;
        5.0 * mt_rate + 3.0 * ox_rate - 5.0 * ens_rate
    }

    fn alpha_like_rate(&self) -> f32 {
        if self.total == 0 {
            return 0.0;
        }
        self.alpha_like_hits as f32 / self.total as f32
    }

    fn unique_rate(&self) -> f32 {
        if self.total == 0 {
            return 0.0;
        }
        self.unique.len() as f32 / self.total as f32
    }
}

/// Detect the gene symbol column in a features.tsv-like file.
pub fn detect_gene_symbol_column(
    reader: &mut dyn Read,
    max_lines: usize,
) -> Result<usize, InputError> {
    let mut reader = BufReader::new(reader);
    let mut stats: Vec<ColumnStats> = Vec::new();
    let mut seen = 0usize;
    let mut line = String::new();

    while seen < max_lines {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|_| InputError::GeneSymbolNotDetected)?;
        if bytes == 0 {
            break;
        }
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.trim_end().split('\t').collect();
        if stats.len() < parts.len() {
            stats.resize_with(parts.len(), ColumnStats::new);
        }
        for (idx, part) in parts.iter().enumerate() {
            let value = part.trim();
            if value.is_empty() {
                continue;
            }
            let col = &mut stats[idx];
            col.total += 1;
            if is_mtdna(value) {
                col.mt_hits += 1;
            }
            if is_oxphos(value) {
                col.oxphos_hits += 1;
            }
            if is_ensembl(value) {
                col.ens_hits += 1;
            }
            if is_alpha_like_symbol(value) {
                col.alpha_like_hits += 1;
            }
            col.unique.insert(value.to_string());
        }
        seen += 1;
    }

    if stats.is_empty() {
        error!("No gene-symbol-like column detected in features.tsv");
        return Err(InputError::GeneSymbolNotDetected);
    }

    let mut best_idx = 0usize;
    let mut best_score = stats[0].score();
    for (idx, stat) in stats.iter().enumerate().skip(1) {
        let score = stat.score();
        if score > best_score {
            best_score = score;
            best_idx = idx;
        }
    }

    if best_score < 0.5 {
        // Fallback for datasets where marker genes were prefiltered:
        // prefer a mostly alphabetic, diverse, non-Ensembl column.
        let fallback = stats
            .iter()
            .enumerate()
            .filter(|(_, s)| s.total > 0)
            .filter(|(_, s)| s.alpha_like_rate() >= 0.7)
            .filter(|(_, s)| s.unique_rate() >= 0.5)
            .filter(|(_, s)| (s.ens_hits as f32 / s.total as f32) < 0.8)
            .max_by(|(_, a), (_, b)| a.unique_rate().total_cmp(&b.unique_rate()))
            .map(|(idx, _)| idx);
        if let Some(idx) = fallback {
            info!(
                index = idx,
                reason = "alphabetic+diverse fallback",
                "Detected gene symbol column"
            );
            return Ok(idx);
        }
        error!(
            "No gene-symbol-like column detected in features.tsv. Hint: file may contain Ensembl IDs only."
        );
        return Err(InputError::GeneSymbolNotDetected);
    }

    let reason = if stats[best_idx].mt_hits > 0 && stats[best_idx].oxphos_hits > 0 {
        "mtDNA+OXPHOS enrichment"
    } else if stats[best_idx].mt_hits > 0 {
        "mtDNA enrichment"
    } else if stats[best_idx].oxphos_hits > 0 {
        "OXPHOS enrichment"
    } else {
        "low marker enrichment"
    };

    info!(
        index = best_idx,
        score = best_score,
        reason,
        "Detected gene symbol column"
    );

    Ok(best_idx)
}

fn is_mtdna(value: &str) -> bool {
    if !value.starts_with("MT-") {
        return false;
    }
    let rest = &value[3..];
    let prefixes = ["ND", "CO", "ATP", "CYB"];
    for prefix in prefixes {
        if rest.starts_with(prefix) {
            let tail = &rest[prefix.len()..];
            return tail.chars().all(|c| c.is_ascii_alphanumeric());
        }
    }
    false
}

fn is_oxphos(value: &str) -> bool {
    value.starts_with("NDUF")
        || value.starts_with("COX")
        || value.starts_with("ATP5")
        || value.starts_with("SDH")
}

fn is_ensembl(value: &str) -> bool {
    if !value.starts_with("ENS") {
        return false;
    }
    let rest = &value[3..];
    let mut idx = 0usize;
    for ch in rest.chars() {
        if ch.is_ascii_uppercase() {
            idx += 1;
            continue;
        }
        break;
    }
    let rest = &rest[idx..];
    if !rest.starts_with('G') {
        return false;
    }
    let digits = &rest[1..];
    !digits.is_empty() && digits.chars().all(|c| c.is_ascii_digit())
}

fn is_alpha_like_symbol(value: &str) -> bool {
    let has_alpha = value.chars().any(|c| c.is_ascii_alphabetic());
    let allowed = value
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.');
    has_alpha && allowed
}

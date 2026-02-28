//! Input format auto-detection.

use std::path::Path;

use kira_scio::detect::DetectedFormat;

use crate::input::{InputError, InputFormat};

/// Concrete input format selected for loading.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectedInputFormat {
    Tenx,
    BDRhapsodyDense,
}

/// Resolve the input format from a user request and path inspection.
pub fn detect_input_format(
    path: &Path,
    requested: InputFormat,
) -> Result<DetectedInputFormat, InputError> {
    match requested {
        InputFormat::Tenx => Ok(DetectedInputFormat::Tenx),
        InputFormat::BdRhapsody => Ok(DetectedInputFormat::BDRhapsodyDense),
        InputFormat::Auto => {
            let detected =
                kira_scio::detect_input_format(path).map_err(|e| InputError::MatrixParse {
                    path: path.to_path_buf(),
                    message: e.message,
                })?;
            Ok(match detected {
                DetectedFormat::BdRhapsodyWta => DetectedInputFormat::BDRhapsodyDense,
                _ => DetectedInputFormat::Tenx,
            })
        }
    }
}

/// Backward-compatible helper kept for existing tests.
pub fn is_probably_bd_rhapsody_dense(path: &Path) -> Result<bool, InputError> {
    let detected = kira_scio::detect_input_format(path).map_err(|e| InputError::MatrixParse {
        path: path.to_path_buf(),
        message: e.message,
    })?;
    Ok(detected == DetectedFormat::BdRhapsodyWta)
}

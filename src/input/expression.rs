//! Canonical expression contract.

use serde::Serialize;

/// Canonical expectations for expression values consumed by kira-mitoqc.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExpressionContract {
    pub normalized: bool,
    pub log1p: bool,
    pub unit: ExpressionUnit,
}

/// Expression unit for normalized values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExpressionUnit {
    CPM,
    TPM,
    Unknown,
}

/// Expression source semantics for downstream stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ExpressionSource {
    RawUmiCounts,
    NormalizedExpression,
}

impl ExpressionContract {
    /// Default v1 contract: log1p CPM/TPM normalized values.
    pub fn v1_default() -> Self {
        Self {
            normalized: true,
            log1p: true,
            unit: ExpressionUnit::Unknown,
        }
    }
}

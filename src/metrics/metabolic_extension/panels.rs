pub const MITO_METABOLIC_PANEL_V1: &str = "MITO_METABOLIC_PANEL_V1";

pub const OXPHOS_PANEL: &[&str] = &[
    "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9",
    "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5",
    "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11", "SDHA", "SDHB", "SDHC", "SDHD",
    "UQCRC1", "UQCRC2", "COX4I1", "COX5A", "COX6C", "ATP5F1A", "ATP5F1B", "ATP5F1C",
];

pub const GLYCOLYSIS_PANEL: &[&str] = &[
    "HK1", "HK2", "PFKP", "PFKM", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA",
];

pub const FAO_PANEL: &[&str] = &["CPT1A", "CPT2", "ACADVL", "ACADM", "HADHA", "HADHB"];

pub const ROS_PANEL: &[&str] = &[
    "SOD1", "SOD2", "GPX1", "GPX4", "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6", "TXN",
    "TXNRD1", "NFE2L2",
];

pub const BIOGENESIS_PANEL: &[&str] =
    &["PPARGC1A", "TFAM", "NRF1", "OPA1", "MFN1", "MFN2", "DNM1L"];

pub fn panel_alias(symbol: &str) -> Option<&'static str> {
    match symbol {
        "NFE2L2" => Some("NRF2"),
        "DNM1L" => Some("DRP1"),
        _ => None,
    }
}

pub fn to_mouse_like(symbol: &str) -> String {
    if symbol.is_empty() {
        return String::new();
    }
    let mut chars = symbol.chars();
    let first = chars.next().map(|c| c.to_ascii_uppercase()).unwrap_or(' ');
    let rest = chars.as_str().to_ascii_lowercase();
    format!("{first}{rest}")
}

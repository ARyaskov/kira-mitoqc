//! Numeric helpers with deterministic behavior.

/// Clamp to [0, 1]. NaN maps to 0.
#[inline]
pub fn clamp01(value: f32) -> f32 {
    if value.is_nan() {
        return 0.0;
    }
    if value < 0.0 {
        0.0
    } else if value > 1.0 {
        1.0
    } else {
        value
    }
}

/// Safe division with an additive epsilon on the denominator.
pub fn safe_div(num: f32, denom: f32, eps: f32) -> f32 {
    num / (denom + eps)
}

/// Approximate equality test.
pub fn approx_eq(a: f32, b: f32, tol: f32) -> bool {
    (a - b).abs() <= tol
}

//! Star-pattern generation and validation helpers.

use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct StarPatternEntry {
    pub x: usize,
    pub y: usize,
    pub magnitude: f64,
    pub use_in_quality: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub median_brightness: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub median_brightest_pixel_part: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct PatternStarSample {
    pub x: usize,
    pub y: usize,
    pub magnitude: f64,
    pub brightest_pixel_part: f64,
}

pub fn median(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() { return None; }
    values.sort_by(|a, b| a.total_cmp(b));
    let n = values.len();
    Some(if n % 2 == 0 { (values[n / 2 - 1] + values[n / 2]) / 2.0 } else { values[n / 2] })
}

pub fn validate_pattern(pattern: &[StarPatternEntry], position_tolerance_px: f64) -> Result<(), String> {
    if !position_tolerance_px.is_finite() || position_tolerance_px <= 0.0 {
        return Err("star_pattern_position_tolerance_px must be greater than zero".into());
    }
    if !pattern.iter().any(|s| s.use_in_quality) {
        return Err("star pattern must contain at least one star with use_in_quality: true".into());
    }
    if pattern.iter().any(|s| !s.magnitude.is_finite() || s.magnitude < 0.0) {
        return Err("star pattern contains invalid magnitude".into());
    }
    Ok(())
}

/// Stable pseudo-random sample selection. Sorting the paths makes the result
/// independent of filesystem enumeration order.
pub fn select_sample_paths(mut paths: Vec<PathBuf>, count: usize, seed: u64) -> Vec<PathBuf> {
    paths.sort();
    if paths.len() <= count { return paths; }
    // Fisher-Yates with a tiny deterministic PRNG; no global RNG is involved.
    let mut state = seed ^ 0x9e3779b97f4a7c15;
    for i in (1..paths.len()).rev() {
        state ^= state << 7; state ^= state >> 9; state ^= state << 8;
        let j = (state as usize) % (i + 1);
        paths.swap(i, j);
    }
    paths.truncate(count);
    paths.sort();
    paths
}

/// Pick three stars near the image centre while retaining a non-degenerate
/// triangle. Brightness is a tie breaker, not a matching criterion.
pub fn recommend_stars(stars: &[PatternStarSample], width: usize, height: usize) -> Vec<usize> {
    if stars.len() < 3 { return Vec::new(); }
    let cx = width as f64 / 2.0; let cy = height as f64 / 2.0;
    let mut ranked: Vec<usize> = (0..stars.len()).collect();
    ranked.sort_by(|&a, &b| {
        let da = ((stars[a].x as f64-cx).powi(2)+(stars[a].y as f64-cy).powi(2)).sqrt();
        let db = ((stars[b].x as f64-cx).powi(2)+(stars[b].y as f64-cy).powi(2)).sqrt();
        da.total_cmp(&db).then_with(|| stars[b].magnitude.total_cmp(&stars[a].magnitude))
    });
    for &a in &ranked { for &b in &ranked { for &c in &ranked {
        if a == b || a == c || b == c { continue; }
        let area = ((stars[b].x as f64-stars[a].x as f64)*(stars[c].y as f64-stars[a].y as f64)
            - (stars[c].x as f64-stars[a].x as f64)*(stars[b].y as f64-stars[a].y as f64)).abs();
        if area > 1.0 { return vec![a, b, c]; }
    }}}
    Vec::new()
}

pub fn default_pattern_path(session: &Path) -> Option<PathBuf> {
    ["stars_pattern.json", "star_pattern.json", "stars.json"].iter()
        .map(|n| session.join(n)).find(|p| p.is_file())
        .or_else(|| ["stars_pattern.json", "star_pattern.json", "stars.json"].iter().map(PathBuf::from).find(|p| p.is_file()))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn median_is_robust_and_even() {
        let mut v = [10.0, 1.0, 2.0, 1000.0];
        assert_eq!(median(&mut v), Some(6.0));
    }
    #[test] fn sampling_is_deterministic_and_bounded() {
        let p: Vec<PathBuf> = (0..30).map(|i| PathBuf::from(format!("{i:02}.fits"))).collect();
        assert_eq!(select_sample_paths(p.clone(), 20, 7), select_sample_paths(p, 20, 7));
    }
    #[test] fn recommendation_has_three_non_collinear_stars() {
        let s = vec![(0,0),(10,0),(0,10),(50,50)].into_iter().map(|(x,y)| PatternStarSample{x,y,magnitude:100.0,brightest_pixel_part:0.2}).collect::<Vec<_>>();
        assert_eq!(recommend_stars(&s, 50, 50).len(), 3);
    }
    #[test] fn pattern_validation_is_explicit() {
        let p = vec![StarPatternEntry{x:1,y:2,magnitude:3.0,use_in_quality:false,median_brightness:None,median_brightest_pixel_part:None}];
        assert!(validate_pattern(&p, 7.0).is_err());
        assert!(validate_pattern(&[], 0.0).is_err());
    }
}

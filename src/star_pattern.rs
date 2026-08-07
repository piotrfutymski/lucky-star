//! Star-pattern generation and validation helpers.

use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct StarPatternEntry {
    pub x: usize,
    pub y: usize,
    pub magnitude: f64,
    pub use_in_quality: bool,
    /// Median integrated star brightness in ADU.
    #[serde(
        default,
        alias = "median_brightness",
        skip_serializing_if = "Option::is_none"
    )]
    pub median_brightness_adu: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub median_brightest_pixel_part: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct PatternStarSample {
    pub x: usize,
    pub y: usize,
    pub magnitude: f64,
    pub magnitude_adu: f64,
    /// Brightest pixel divided by the detector range (0..1), not a fraction
    /// of the star flux.  Keeping this name preserves the original API.
    pub brightest_pixel_part: f64,
}

#[derive(Clone, Debug)]
pub struct AggregatedStar {
    pub reference_index: usize,
    pub sample: PatternStarSample,
    pub sample_count: usize,
}

pub fn median(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let n = values.len();
    Some(if n.is_multiple_of(2) {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    } else {
        values[n / 2]
    })
}

pub fn validate_pattern(
    pattern: &[StarPatternEntry],
    position_tolerance_px: f64,
) -> Result<(), String> {
    if !position_tolerance_px.is_finite() || position_tolerance_px <= 0.0 {
        return Err("star_pattern_position_tolerance_px must be greater than zero".into());
    }
    if !pattern.iter().any(|s| s.use_in_quality) {
        return Err("star pattern must contain at least one star with use_in_quality: true".into());
    }
    if pattern
        .iter()
        .any(|s| !s.magnitude.is_finite() || s.magnitude < 0.0)
    {
        return Err("star pattern contains invalid magnitude".into());
    }
    Ok(())
}

/// Stable pseudo-random sample selection. Sorting the paths makes the result
/// independent of filesystem enumeration order.
pub fn select_sample_paths(mut paths: Vec<PathBuf>, count: usize, seed: u64) -> Vec<PathBuf> {
    paths.sort();
    if paths.len() <= count {
        return paths;
    }
    let mut state = seed ^ 0x9e3779b97f4a7c15;
    for i in (1..paths.len()).rev() {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        let j = (state as usize) % (i + 1);
        paths.swap(i, j);
    }
    paths.truncate(count);
    paths.sort();
    paths
}

fn transform_point(p: &PatternStarSample, tx: f64, ty: f64, angle: f64) -> (f64, f64) {
    let (sin, cos) = angle.sin_cos();
    (
        cos * p.x as f64 - sin * p.y as f64 + tx,
        sin * p.x as f64 + cos * p.y as f64 + ty,
    )
}

fn distance(a: &PatternStarSample, b: &PatternStarSample) -> f64 {
    ((a.x as f64 - b.x as f64).powi(2) + (a.y as f64 - b.y as f64).powi(2)).sqrt()
}

/// Find a rigid transform using only pairwise geometry.  No brightness or
/// detector-order cutoff is applied: every detected point remains a candidate.
fn best_transform(
    reference: &[PatternStarSample],
    current: &[PatternStarSample],
    tolerance: f64,
) -> Option<(f64, f64, f64, usize, f64)> {
    if reference.len() < 2 || current.len() < 2 {
        return None;
    }
    let mut best: Option<(f64, f64, f64, usize, f64)> = None;
    for a in 0..reference.len() {
        for b in 0..reference.len() {
            if a == b {
                continue;
            }
            let ref_len = distance(&reference[a], &reference[b]);
            if ref_len < 2.0 {
                continue;
            }
            for c in 0..current.len() {
                for d in 0..current.len() {
                    if c == d {
                        continue;
                    }
                    let cur_len = distance(&current[c], &current[d]);
                    if (cur_len - ref_len).abs() > tolerance * 2.5 {
                        continue;
                    }
                    let angle = (current[d].y as f64 - current[c].y as f64)
                        .atan2(current[d].x as f64 - current[c].x as f64)
                        - (reference[b].y as f64 - reference[a].y as f64)
                            .atan2(reference[b].x as f64 - reference[a].x as f64);
                    let (rx, ry) = transform_point(&reference[a], 0.0, 0.0, angle);
                    let tx = current[c].x as f64 - rx;
                    let ty = current[c].y as f64 - ry;
                    let mut count = 0;
                    let mut error = 0.0;
                    for r in reference.iter() {
                        let (x, y) = transform_point(r, tx, ty, angle);
                        let nearest = current
                            .iter()
                            .map(|q| ((q.x as f64 - x).powi(2) + (q.y as f64 - y).powi(2)).sqrt())
                            .fold(f64::INFINITY, f64::min);
                        if nearest <= tolerance {
                            count += 1;
                            error += nearest;
                        }
                    }
                    let candidate = (tx, ty, angle, count, error);
                    if best
                        .map(|b| count > b.3 || (count == b.3 && error < b.4))
                        .unwrap_or(true)
                    {
                        best = Some(candidate);
                    }
                }
            }
        }
    }
    best.filter(|b| b.3 >= 2)
}

/// Aggregate detections from multiple frames into stars in the first frame's
/// coordinate system.  A star is retained only when it is seen in at least
/// half of the sampled frames, and each brightness field is a median.
pub fn aggregate_samples(
    reference: &[PatternStarSample],
    frames: &[Vec<PatternStarSample>],
    tolerance: f64,
) -> Vec<AggregatedStar> {
    if reference.is_empty() || frames.is_empty() {
        return Vec::new();
    }
    let mut tracks: Vec<Vec<PatternStarSample>> =
        reference.iter().cloned().map(|s| vec![s]).collect();
    for frame in frames.iter().skip(1) {
        let Some((tx, ty, angle, _, _)) = best_transform(reference, frame, tolerance) else {
            continue;
        };
        let mut used = vec![false; frame.len()];
        for (ri, ref_star) in reference.iter().enumerate() {
            let (x, y) = transform_point(ref_star, tx, ty, angle);
            let nearest = frame
                .iter()
                .enumerate()
                .filter(|(i, _)| !used[*i])
                .map(|(i, s)| {
                    (
                        i,
                        ((s.x as f64 - x).powi(2) + (s.y as f64 - y).powi(2)).sqrt(),
                    )
                })
                .min_by(|a, b| a.1.total_cmp(&b.1));
            if let Some((i, d)) = nearest
                && d <= tolerance
            {
                used[i] = true;
                tracks[ri].push(frame[i].clone());
            }
        }
    }
    let required = (frames.len().div_ceil(2)).max(1);
    tracks
        .into_iter()
        .enumerate()
        .filter_map(|(reference_index, values)| {
            if values.len() < required {
                return None;
            }
            let mut magnitudes: Vec<f64> = values.iter().map(|s| s.magnitude).collect();
            let mut brightness: Vec<f64> = values.iter().map(|s| s.brightest_pixel_part).collect();
            Some(AggregatedStar {
                reference_index,
                sample: PatternStarSample {
                    x: reference[reference_index].x,
                    y: reference[reference_index].y,
                    magnitude: median(&mut magnitudes)?,
                    magnitude_adu: {
                        let mut brightness_adu: Vec<f64> =
                            values.iter().map(|s| s.magnitude_adu).collect();
                        median(&mut brightness_adu)?
                    },
                    brightest_pixel_part: median(&mut brightness)?,
                },
                sample_count: values.len(),
            })
        })
        .collect()
}

/// Pick three stars near the image centre while retaining a non-degenerate
/// triangle.  Stars in the preferred 10..=50% detector range are ranked
/// before fallback candidates outside that range.
pub fn recommend_stars(stars: &[PatternStarSample], width: usize, height: usize) -> Vec<usize> {
    if stars.len() < 3 {
        return Vec::new();
    }
    let cx = width as f64 / 2.0;
    let cy = height as f64 / 2.0;
    let mut ranked: Vec<usize> = (0..stars.len()).collect();
    ranked.sort_by(|&a, &b| {
        let preferred_a = (0.10..=0.50).contains(&stars[a].brightest_pixel_part);
        let preferred_b = (0.10..=0.50).contains(&stars[b].brightest_pixel_part);
        preferred_b
            .cmp(&preferred_a)
            .then_with(|| {
                let da = (stars[a].x as f64 - cx).powi(2) + (stars[a].y as f64 - cy).powi(2);
                let db = (stars[b].x as f64 - cx).powi(2) + (stars[b].y as f64 - cy).powi(2);
                da.total_cmp(&db)
            })
            .then_with(|| stars[b].magnitude.total_cmp(&stars[a].magnitude))
    });
    let preferred = |i: usize| (0.10..=0.50).contains(&stars[i].brightest_pixel_part);
    let mut best: Option<(usize, f64, Vec<usize>)> = None;
    for &a in &ranked {
        for &b in &ranked {
            for &c in &ranked {
                if a == b || a == c || b == c {
                    continue;
                }
                let area = ((stars[b].x as f64 - stars[a].x as f64)
                    * (stars[c].y as f64 - stars[a].y as f64)
                    - (stars[c].x as f64 - stars[a].x as f64)
                        * (stars[b].y as f64 - stars[a].y as f64))
                    .abs();
                if area <= 1.0 {
                    continue;
                }
                let count = usize::from(preferred(a))
                    + usize::from(preferred(b))
                    + usize::from(preferred(c));
                let centre_distance = [a, b, c]
                    .iter()
                    .map(|i| (stars[*i].x as f64 - cx).powi(2) + (stars[*i].y as f64 - cy).powi(2))
                    .sum::<f64>();
                if best
                    .as_ref()
                    .map(|(old_count, old_distance, _)| {
                        count > *old_count
                            || (count == *old_count && centre_distance < *old_distance)
                    })
                    .unwrap_or(true)
                {
                    best = Some((count, centre_distance, vec![a, b, c]));
                }
            }
        }
    }
    best.map(|(_, _, indices)| indices).unwrap_or_default()
}

pub fn default_pattern_path(session: &Path) -> Option<PathBuf> {
    ["stars_pattern.json", "star_pattern.json", "stars.json"]
        .iter()
        .map(|n| session.join(n))
        .find(|p| p.is_file())
        .or_else(|| {
            ["stars_pattern.json", "star_pattern.json", "stars.json"]
                .iter()
                .map(PathBuf::from)
                .find(|p| p.is_file())
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn sample(x: usize, y: usize, m: f64, b: f64) -> PatternStarSample {
        PatternStarSample {
            x,
            y,
            magnitude: m,
            magnitude_adu: m * 10.0,
            brightest_pixel_part: b,
        }
    }
    #[test]
    fn median_is_robust_and_even() {
        let mut v = [10.0, 1.0, 2.0, 1000.0];
        assert_eq!(median(&mut v), Some(6.0));
    }
    #[test]
    fn sampling_is_deterministic_and_bounded() {
        let p: Vec<PathBuf> = (0..30)
            .map(|i| PathBuf::from(format!("{i:02}.fits")))
            .collect();
        assert_eq!(
            select_sample_paths(p.clone(), 20, 7),
            select_sample_paths(p, 20, 7)
        );
    }
    #[test]
    fn recommendation_prefers_detector_range() {
        let s = vec![
            sample(5, 5, 100.0, 0.01),
            sample(45, 5, 90.0, 0.2),
            sample(5, 45, 80.0, 0.3),
            sample(45, 45, 70.0, 0.4),
        ];
        let result = recommend_stars(&s, 50, 50);
        assert!(result.iter().all(|i| *i == 1 || *i == 2 || *i == 3));
        assert_eq!(result.len(), 3);
    }
    #[test]
    fn aggregation_uses_medians_and_reference_numbers() {
        let reference = vec![
            sample(10, 10, 10.0, 0.2),
            sample(30, 10, 20.0, 0.2),
            sample(10, 30, 30.0, 0.2),
        ];
        let frames = vec![
            reference.clone(),
            vec![
                sample(11, 10, 100.0, 0.4),
                sample(31, 10, 200.0, 0.4),
                sample(11, 30, 300.0, 0.4),
            ],
        ];
        let a = aggregate_samples(&reference, &frames, 3.0);
        assert_eq!(a.len(), 3);
        assert_eq!(a[0].sample.magnitude, 55.0);
        assert_eq!(a[0].reference_index, 0);
    }
    #[test]
    fn pattern_validation_is_explicit() {
        let p = vec![StarPatternEntry {
            x: 1,
            y: 2,
            magnitude: 3.0,
            use_in_quality: false,
            median_brightness_adu: None,
            median_brightest_pixel_part: None,
        }];
        assert!(validate_pattern(&p, 7.0).is_err());
        assert!(validate_pattern(&[], 0.0).is_err());
    }
}

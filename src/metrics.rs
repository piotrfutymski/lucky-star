//! Independent metric values and filtering primitives.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Clone, Debug, PartialEq)]
pub struct MetricValues {
    pub quality: f64, pub fwhm: f64, pub star_count: usize,
    pub brightest_star_adu: f64, pub brightest_star_photons: f64, pub star5_photons: f64,
    pub background_raw_adu: f64, pub background_corrected_adu: f64,
    pub quality_star_pattern: Option<f64>, pub quality_star_pattern_source: bool,
    pub star_brightness_adu: Option<f64>, pub snr: Option<f64>,
    pub star_pattern_found: bool, pub matched_star_count: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MetricRecord { pub file_name: String, pub modified_ns: u128, pub metrics: MetricValues }

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Metric {
    Quality,
    Fwhm,
    QualityStarPattern,
    Background,
    StarBrightness,
    Snr,
}

impl Metric {
    pub fn key(self) -> &'static str {
        match self {
            Self::Quality => "quality",
            Self::Fwhm => "fwhm",
            Self::QualityStarPattern => "quality_star_pattern",
            Self::Background => "background",
            Self::StarBrightness => "star_brightness",
            Self::Snr => "snr",
        }
    }

    pub fn higher_is_better(self) -> bool {
        !matches!(self, Self::Background | Self::Fwhm)
    }

    pub fn requires_pattern(self) -> bool {
        matches!(
            self,
            Self::QualityStarPattern | Self::StarBrightness | Self::Snr
        )
    }

    pub fn value(self, m: &MetricValues) -> Option<f64> {
        match self {
            Self::Quality => Some(m.quality),
            Self::Fwhm => Some(m.fwhm),
            Self::QualityStarPattern => m.quality_star_pattern.or(Some(m.quality)),
            Self::Background => Some(m.background_corrected_adu),
            Self::StarBrightness => m.star_brightness_adu,
            Self::Snr => m.snr,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FilterRule {
    pub metric: Metric,
    pub relative: Option<f64>,
    pub absolute: Option<f64>,
}

impl FilterRule {
    pub fn relative(metric: Metric, multiplier: f64) -> Result<Self, String> {
        if !multiplier.is_finite() || multiplier <= 0.0 {
            return Err("relative threshold must be greater than zero".into());
        }
        Ok(Self {
            metric,
            relative: Some(multiplier),
            absolute: None,
        })
    }

    pub fn absolute(metric: Metric, threshold: f64) -> Result<Self, String> {
        if !threshold.is_finite() {
            return Err("absolute threshold must be finite".into());
        }
        Ok(Self {
            metric,
            relative: None,
            absolute: Some(threshold),
        })
    }

    pub fn threshold(&self, median: Option<f64>) -> Result<f64, String> {
        match (self.relative, self.absolute) {
            (Some(_), Some(_)) => Err(format!(
                "both relative and absolute thresholds supplied for {}",
                self.metric.key()
            )),
            (Some(factor), None) => median
                .map(|m| factor * m)
                .ok_or_else(|| format!("no values available for {}", self.metric.key())),
            (None, Some(value)) => Ok(value),
            (None, None) => Err("filter rule has no threshold".into()),
        }
    }
}

pub fn median(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    Some(if values.len().is_multiple_of(2) {
        (values[values.len() / 2 - 1] + values[values.len() / 2]) / 2.0
    } else {
        values[values.len() / 2]
    })
}

pub fn medians(records: &[MetricRecord]) -> BTreeMap<Metric, f64> {
    let mut result = BTreeMap::new();
    for metric in [
        Metric::Quality,
        Metric::Fwhm,
        Metric::QualityStarPattern,
        Metric::Background,
        Metric::StarBrightness,
        Metric::Snr,
    ] {
        let mut values: Vec<f64> = records
            .iter()
            .filter_map(|r| if metric == Metric::QualityStarPattern && !r.metrics.quality_star_pattern_source { None } else { metric.value(&r.metrics) })
            .filter(|v| v.is_finite())
            .collect();
        if let Some(m) = median(&mut values) {
            result.insert(metric, m);
        }
    }
    result
}

pub fn pattern_quality_medians(records: &[MetricRecord]) -> (Option<f64>, Option<f64>) {
    let mut primary: Vec<f64> = records.iter().filter(|r| r.metrics.quality_star_pattern_source)
        .filter_map(|r| r.metrics.quality_star_pattern).filter(|v| v.is_finite()).collect();
    let mut fallback: Vec<f64> = records.iter().filter(|r| !r.metrics.quality_star_pattern_source)
        .map(|r| r.metrics.quality).filter(|v| v.is_finite()).collect();
    (median(&mut primary), median(&mut fallback))
}

pub fn passes_all_filters(
    record: &MetricRecord,
    rules: &[FilterRule],
    medians: &BTreeMap<Metric, f64>,
) -> Result<bool, String> {
    for rule in rules {
        let Some(value) = rule.metric.value(&record.metrics) else {
            return Ok(false);
        };
        let median_value = if rule.metric == Metric::QualityStarPattern && rule.relative.is_some() && !record.metrics.quality_star_pattern_source {
            medians.get(&Metric::Quality).copied()
        } else { medians.get(&rule.metric).copied() };
        let threshold = rule.threshold(median_value)?;
        let passes = if rule.metric.higher_is_better() {
            value >= threshold
        } else {
            value <= threshold
        };
        if !passes {
            return Ok(false);
        }
    }
    Ok(true)
}

pub fn unique_removed_folder(dir: &Path, rules: &[FilterRule], timestamp: u64) -> PathBuf {
    let base = if rules.len() == 1 {
        let r = &rules[0];
        let suffix = r
            .relative
            .map(|v| format!("{v:.3}"))
            .or_else(|| r.absolute.map(|v| format!("{v:.3}")))
            .unwrap_or_default();
        dir.join(format!("removed_{}_{}", r.metric.key(), suffix))
    } else {
        dir.join(format!("removed_combined_{timestamp}"))
    };
    if !base.exists() {
        return base;
    }
    for n in 2.. {
        let candidate = PathBuf::from(format!("{}_{}", base.display(), n));
        if !candidate.exists() {
            return candidate;
        }
    }
    unreachable!()
}

pub fn current_timestamp() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs()
}


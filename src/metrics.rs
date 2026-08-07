//! Phase 2: durable metrics, cache invalidation and filtering primitives.
//!
//! The FITS-specific adapter lives in `main.rs`; this module deliberately keeps
//! the cache and filtering rules independent of the command line and renderer.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

pub const CACHE_FORMAT_VERSION: u32 = 1;
pub const ALGORITHM_VERSION: &str = "phase2-metrics-1";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct MetricValues {
    pub quality: f64,
    pub fwhm: f64,
    pub star_count: usize,
    pub brightest_star_adu: f64,
    pub background_raw_adu: f64,
    pub background_corrected_adu: f64,
    pub quality_star_pattern: Option<f64>,
    pub star_brightness_adu: Option<f64>,
    pub snr: Option<f64>,
    pub star_pattern_found: bool,
    pub matched_star_count: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct CacheRecord {
    pub file_name: String,
    pub size: u64,
    pub modified_ns: u128,
    pub cache_format_version: u32,
    pub algorithm_version: String,
    pub configuration_fingerprint: String,
    pub metrics: MetricValues,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct MetricsCache {
    pub cache_format_version: u32,
    pub algorithm_version: String,
    pub configuration_fingerprint: String,
    pub records: Vec<CacheRecord>,
}

impl MetricsCache {
    pub fn empty(fingerprint: impl Into<String>) -> Self {
        Self {
            cache_format_version: CACHE_FORMAT_VERSION,
            algorithm_version: ALGORITHM_VERSION.into(),
            configuration_fingerprint: fingerprint.into(),
            records: Vec::new(),
        }
    }

    pub fn is_compatible(&self, fingerprint: &str) -> bool {
        self.cache_format_version == CACHE_FORMAT_VERSION
            && self.algorithm_version == ALGORITHM_VERSION
            && self.configuration_fingerprint == fingerprint
            && self.records.iter().all(|r| {
                r.cache_format_version == CACHE_FORMAT_VERSION
                    && r.algorithm_version == ALGORITHM_VERSION
                    && r.configuration_fingerprint == fingerprint
            })
    }

    pub fn load(path: &Path, fingerprint: &str) -> Result<Self, String> {
        let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
        let cache: Self = serde_json::from_str(&text).map_err(|e| e.to_string())?;
        Ok(if cache.is_compatible(fingerprint) {
            cache
        } else {
            Self::empty(fingerprint)
        })
    }

    pub fn save(&self, path: &Path) -> Result<(), String> {
        let text = serde_json::to_string_pretty(self).map_err(|e| e.to_string())?;
        fs::write(path, text).map_err(|e| e.to_string())
    }

    pub fn upsert(&mut self, record: CacheRecord) {
        if let Some(old) = self
            .records
            .iter_mut()
            .find(|r| r.file_name == record.file_name)
        {
            *old = record;
        } else {
            self.records.push(record);
        }
        self.records.sort_by(|a, b| a.file_name.cmp(&b.file_name));
    }

    pub fn record_is_current(&self, path: &Path, fingerprint: &str) -> bool {
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            return false;
        };
        let Ok(meta) = fs::metadata(path) else {
            return false;
        };
        self.records.iter().any(|r| {
            r.file_name == name
                && r.size == meta.len()
                && r.modified_ns == modified_ns(&meta)
                && r.configuration_fingerprint == fingerprint
                && r.algorithm_version == ALGORITHM_VERSION
                && r.cache_format_version == CACHE_FORMAT_VERSION
        })
    }
}

pub fn modified_ns(meta: &fs::Metadata) -> u128 {
    meta.modified()
        .ok()
        .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
        .map(|d| d.as_nanos())
        .unwrap_or(0)
}

pub fn record_for(
    path: &Path,
    fingerprint: impl Into<String>,
    metrics: MetricValues,
) -> Result<CacheRecord, String> {
    let meta = fs::metadata(path).map_err(|e| e.to_string())?;
    Ok(CacheRecord {
        file_name: path
            .file_name()
            .and_then(|n| n.to_str())
            .ok_or("invalid file name")?
            .to_string(),
        size: meta.len(),
        modified_ns: modified_ns(&meta),
        cache_format_version: CACHE_FORMAT_VERSION,
        algorithm_version: ALGORITHM_VERSION.into(),
        configuration_fingerprint: fingerprint.into(),
        metrics,
    })
}

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
            Self::QualityStarPattern => m.quality_star_pattern,
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
    Some(if values.len() % 2 == 0 {
        (values[values.len() / 2 - 1] + values[values.len() / 2]) / 2.0
    } else {
        values[values.len() / 2]
    })
}

pub fn medians(records: &[CacheRecord]) -> BTreeMap<Metric, f64> {
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
            .filter_map(|r| metric.value(&r.metrics))
            .filter(|v| v.is_finite())
            .collect();
        if let Some(m) = median(&mut values) {
            result.insert(metric, m);
        }
    }
    result
}

pub fn passes_all_filters(
    record: &CacheRecord,
    rules: &[FilterRule],
    medians: &BTreeMap<Metric, f64>,
) -> Result<bool, String> {
    for rule in rules {
        let Some(value) = rule.metric.value(&record.metrics) else {
            return Ok(false);
        };
        let threshold = rule.threshold(medians.get(&rule.metric).copied())?;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn values(quality: f64, snr: Option<f64>) -> MetricValues {
        MetricValues {
            quality,
            fwhm: 2.0,
            star_count: 3,
            brightest_star_adu: 10.0,
            background_raw_adu: 20.0,
            background_corrected_adu: 20.0,
            quality_star_pattern: snr,
            star_brightness_adu: snr.map(|x| x * 10.0),
            snr,
            star_pattern_found: snr.is_some(),
            matched_star_count: usize::from(snr.is_some()),
        }
    }

    fn record(name: &str, quality: f64, snr: Option<f64>) -> CacheRecord {
        CacheRecord {
            file_name: name.into(),
            size: 1,
            modified_ns: 1,
            cache_format_version: CACHE_FORMAT_VERSION,
            algorithm_version: ALGORITHM_VERSION.into(),
            configuration_fingerprint: "config-a".into(),
            metrics: values(quality, snr),
        }
    }

    #[test]
    fn cache_round_trip_and_stale_configuration_invalidates_everything() {
        let dir = std::env::temp_dir().join(format!("lucky-star-phase2-{}", current_timestamp()));
        let _ = fs::create_dir_all(&dir);
        let path = dir.join("metrics_cache.json");
        let mut cache = MetricsCache::empty("config-a");
        cache.upsert(record("a.fits", 0.9, Some(10.0)));
        cache.save(&path).unwrap();
        assert_eq!(
            MetricsCache::load(&path, "config-a").unwrap().records.len(),
            1
        );
        assert!(
            MetricsCache::load(&path, "config-b")
                .unwrap()
                .records
                .is_empty()
        );
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn relative_rules_use_median_and_are_and_combined() {
        let records = vec![
            record("a", 1.0, Some(10.0)),
            record("b", 0.4, Some(20.0)),
            record("c", 0.6, None),
        ];
        let med = medians(&records);
        let quality = FilterRule::relative(Metric::Quality, 0.83).unwrap();
        let snr = FilterRule::relative(Metric::Snr, 0.707).unwrap();
        assert!(passes_all_filters(&records[0], &[quality.clone()], &med).unwrap());
        assert!(!passes_all_filters(&records[1], &[quality, snr.clone()], &med).unwrap());
        assert!(!passes_all_filters(&records[2], &[snr], &med).unwrap());
    }

    #[test]
    fn absolute_rules_obey_metric_direction_and_conflicts_are_rejected() {
        let r = record("a", 0.8, Some(10.0));
        let med = medians(std::slice::from_ref(&r));
        assert!(
            passes_all_filters(
                &r,
                &[FilterRule::absolute(Metric::Fwhm, 2.0).unwrap()],
                &med
            )
            .unwrap()
        );
        assert!(
            !passes_all_filters(
                &r,
                &[FilterRule::absolute(Metric::Background, 1.0).unwrap()],
                &med
            )
            .unwrap()
        );
        let mut both = FilterRule::relative(Metric::Snr, 1.0).unwrap();
        both.absolute = Some(2.0);
        assert!(both.threshold(Some(10.0)).is_err());
    }

    #[test]
    fn cache_records_keep_missing_pattern_metrics_as_missing() {
        let r = record("cloudy.fits", 0.5, None);
        assert_eq!(r.metrics.quality_star_pattern, None);
        assert_eq!(r.metrics.star_brightness_adu, None);
        assert_eq!(r.metrics.snr, None);
        assert!(!r.metrics.star_pattern_found);
        assert!(
            !passes_all_filters(
                &r,
                &[FilterRule::absolute(Metric::Snr, 1.0).unwrap()],
                &BTreeMap::new()
            )
            .unwrap()
        );
    }

    #[test]
    fn removed_folder_never_overwrites_an_existing_folder() {
        let dir = std::env::temp_dir().join(format!("lucky-star-folder-{}", current_timestamp()));
        fs::create_dir_all(dir.join("removed_snr_0.707")).unwrap();
        let rule = FilterRule::relative(Metric::Snr, 0.707).unwrap();
        let result = unique_removed_folder(&dir, &[rule], 42);
        assert_ne!(result, dir.join("removed_snr_0.707"));
        let _ = fs::remove_dir_all(dir);
    }
}

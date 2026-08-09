use crate::astro_image::{AstroImage, fwhm_from_quality};
use crate::constellation::{Constellation, RegisteredStar, load_stars_from_json};
use crate::star::Star;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

pub mod astro_image;
pub mod constellation;
pub mod helpers;
pub mod metrics;
pub mod star;
pub mod star_pattern;

#[derive(Serialize)]
pub struct AppConfig {
    gain_to_adu: HashMap<u32, f64>,
    min_photons_to_detect_star: i32,
    min_central_photons_to_detect_star: i32,
    psf_size: usize,
    min_photons_quality: f64,
    rolling_avg_window: usize,
    log_quality_window_t: f64,
    pub star_pattern_position_tolerance_px: f64,
    pub background_bias_adu: f64,
}

impl Default for AppConfig {
    fn default() -> Self {
        let mut gain_to_adu = HashMap::new();
        gain_to_adu.insert(5200, 27.5);
        gain_to_adu.insert(7000, 40.0);
        gain_to_adu.insert(9000, 61.5);
        gain_to_adu.insert(12000, 143.0);
        gain_to_adu.insert(15000, 250.0);
        gain_to_adu.insert(0, 250.0);
        AppConfig {
            gain_to_adu,
            min_photons_to_detect_star: 300,
            min_central_photons_to_detect_star: 30,
            psf_size: 13,
            min_photons_quality: 450.0,
            rolling_avg_window: 100,
            log_quality_window_t: 300.0,
            star_pattern_position_tolerance_px: 7.0,
            background_bias_adu: 1000.0,
        }
    }
}

impl<'de> Deserialize<'de> for AppConfig {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct PartialAppConfig {
            #[serde(default)]
            gain_to_adu: Option<HashMap<u32, f64>>,
            #[serde(default)]
            min_photons_to_detect_star: Option<i32>,
            #[serde(default)]
            min_central_photons_to_detect_star: Option<i32>,
            #[serde(default)]
            psf_size: Option<usize>,
            #[serde(default)]
            min_photons_quality: Option<f64>,
            #[serde(default)]
            rolling_avg_window: Option<usize>,
            #[serde(default)]
            log_quality_window_t: Option<f64>,
            #[serde(default)]
            star_pattern_position_tolerance_px: Option<f64>,
            background_bias_adu: Option<f64>,
        }

        let partial = PartialAppConfig::deserialize(deserializer)?;
        let default = AppConfig::default();

        Ok(AppConfig {
            gain_to_adu: partial.gain_to_adu.unwrap_or(default.gain_to_adu),
            min_photons_to_detect_star: partial
                .min_photons_to_detect_star
                .unwrap_or(default.min_photons_to_detect_star),
            min_central_photons_to_detect_star: partial
                .min_central_photons_to_detect_star
                .unwrap_or(default.min_central_photons_to_detect_star),
            psf_size: partial.psf_size.unwrap_or(default.psf_size),
            min_photons_quality: partial
                .min_photons_quality
                .unwrap_or(default.min_photons_quality),
            rolling_avg_window: partial
                .rolling_avg_window
                .unwrap_or(default.rolling_avg_window),
            log_quality_window_t: partial
                .log_quality_window_t
                .unwrap_or(default.log_quality_window_t),
            star_pattern_position_tolerance_px: partial
                .star_pattern_position_tolerance_px
                .unwrap_or(default.star_pattern_position_tolerance_px),
            background_bias_adu: partial
                .background_bias_adu
                .unwrap_or(default.background_bias_adu),
        })
    }
}

impl AppConfig {
    fn validate(&self) -> Result<(), String> {
        if !self.background_bias_adu.is_finite() || self.background_bias_adu < 0.0 {
            return Err("background_bias_adu must be finite and must not be negative".into());
        }
        if !self.star_pattern_position_tolerance_px.is_finite()
            || self.star_pattern_position_tolerance_px <= 0.0
        {
            return Err("star_pattern_position_tolerance_px must be greater than zero".into());
        }
        Ok(())
    }

    fn load_or_default() -> Self {
        // Try executable directory first
        if let Ok(exe_path) = std::env::current_exe()
            && let Some(exe_dir) = exe_path.parent()
        {
            let config_path = exe_dir.join("config.json");
            if let Ok(content) = fs::read_to_string(&config_path) {
                match serde_json::from_str(&content) {
                    Ok(cfg) => return cfg,
                    Err(e) => eprintln!(
                        "Warning: failed to parse config.json in executable directory: {}",
                        e
                    ),
                }
            }
        }
        // Try current working directory
        if let Ok(cwd) = std::env::current_dir() {
            let config_path = cwd.join("config.json");
            if let Ok(content) = fs::read_to_string(&config_path) {
                match serde_json::from_str(&content) {
                    Ok(cfg) => return cfg,
                    Err(e) => eprintln!(
                        "Warning: failed to parse config.json in current directory: {}",
                        e
                    ),
                }
            }
        }
        AppConfig::default()
    }
}

#[derive(Parser)]
#[command(name = "lucky-star", about = "Astronomical image quality analyzer")]
struct Args {
    /// Path to a FITS file or a directory containing FITS files
    #[arg(default_value = ".")]
    path: String,

    /// Apply metric filters and move rejected FITS files to a new folder.
    #[arg(long)]
    filter: bool,

    /// After filtering, divide kept files into batches: --divide METRIC COUNT
    #[arg(long, value_names = ["METRIC", "COUNT"], num_args = 2)]
    divide: Option<Vec<String>>,

    /// Show detailed processing diagnostics.
    #[arg(long)]
    verbose: bool,

    #[arg(long)]
    snr: Option<f64>,
    #[arg(long)]
    quality_star_pattern: Option<f64>,
    #[arg(long)]
    background: Option<f64>,
    #[arg(long)]
    star_brightness: Option<f64>,
    #[arg(long)]
    quality: Option<f64>,
    #[arg(long)]
    fwhm: Option<f64>,
    #[arg(long)]
    snr_absolute: Option<f64>,
    #[arg(long)]
    background_absolute: Option<f64>,
    #[arg(long)]
    star_brightness_absolute: Option<f64>,
    #[arg(long)]
    quality_star_pattern_absolute: Option<f64>,
    #[arg(long)]
    quality_absolute: Option<f64>,
    #[arg(long)]
    fwhm_absolute: Option<f64>,

    /// Only search for stars in the central fraction of the image (e.g. 0.3 = central 30% width and height)
    #[arg(long, value_name = "FRACTION")]
    crop: Option<f64>,

    /// Save annotated star image to a JPG file
    #[arg(long, short)]
    save_stars: bool,

    /// Path to a JSON file with reference stars for constellation-based quality filtering
    #[arg(long, value_name = "FILE")]
    star_pattern: Option<String>,

    /// Quick analysis of the N newest FITS images; skips charts.
    #[arg(short = 'c', long = "check-count", value_name = "N")]
    check_count: Option<usize>,

    /// Interactively generate stars_pattern.json from a session folder
    #[arg(long, value_name = "FOLDER")]
    make_star_pattern: Option<String>,

    /// Comma-separated star numbers from star_pattern_candidates.jpg
    #[arg(long, value_name = "N,N,N")]
    star_pattern_numbers: Option<String>,
}

struct ImageInfo {
    file_name: String,
    file_path: PathBuf,
    quality: f64,
    fwhm: f64,
    quality_image: Option<f64>,
    star_count: usize,
    constellation_found: Option<bool>,
    matched_star_count: usize,
    star_brightness_adu: Option<f64>,
    snr: Option<f64>,
    brightest_star_adu: f64,
    background_raw_adu: f64,
    background_corrected_adu: f64,
    brightest_star_photons: f64,
    star5_photons: f64,
}

/// Sum only the registered stars explicitly enabled for quality metrics.
/// This is also the signal used to derive the star-pattern SNR.
fn star_pattern_brightness_adu(stars: &[Star], constellation: &Constellation) -> Option<f64> {
    if !constellation.found {
        return None;
    }
    Some(
        constellation
            .registered_stars
            .iter()
            .enumerate()
            .filter(|(_, star)| star.use_in_quality)
            .filter_map(|(registered_index, _)| {
                constellation
                    .star_mapping
                    .get(&registered_index)
                    .and_then(|detected_index| stars.get(*detected_index))
            })
            .map(|star| star.magnitude_adu)
            .sum(),
    )
}

fn star_pattern_snr(signal: Option<f64>, background_noise_adu: f64) -> Option<f64> {
    signal.and_then(|signal| {
        (background_noise_adu.is_finite() && background_noise_adu > 0.0)
            .then_some(signal / background_noise_adu)
    })
}

fn apply_constellation_quality(
    img: &mut AstroImage,
    registered_stars: &[RegisteredStar],
    config: &AppConfig,
    label: &str,
) -> bool {
    let constellation = Constellation::find_in_image_with_tolerance(
        registered_stars.to_vec(),
        img,
        config.star_pattern_position_tolerance_px as f32,
    );
    if constellation.found {
        let quality_indices: HashSet<usize> = constellation
            .registered_stars
            .iter()
            .enumerate()
            .filter(|(_, rs)| rs.use_in_quality)
            .filter_map(|(i, _)| constellation.star_mapping.get(&i).copied())
            .collect();
        img.recalculate_quality_for_star_indices(&quality_indices, config);
        true
    } else {
        eprintln!("Warning: star-pattern metrics unavailable for '{}'.", label);
        false
    }
}

fn process_single_file(
    path: &Path,
    crop: Option<f64>,
    save_stars: bool,
    config: &AppConfig,
    registered_stars: Option<&[RegisteredStar]>,
) {
    let mut img = AstroImage::load(path, crop, config).unwrap();
    if let Some(stars) = registered_stars {
        apply_constellation_quality(&mut img, stars, config, &path.display().to_string());
    }
    print!("{}", img);
    if save_stars {
        let stem = path.file_stem().unwrap_or_default().to_string_lossy();
        let dir = path.parent().unwrap_or_else(|| Path::new("."));
        let jpg_path = dir.join(format!("{}_stars.jpg", stem));
        let md_path = dir.join(format!("{}_stars.md", stem));
        img.save_stars_jpg(&jpg_path).unwrap();
        img.save_stars_md(&md_path).unwrap();
        println!("Stars image saved to: {}", jpg_path.display());
        println!("Stars table saved to: {}", md_path.display());
    }
}

fn collect_fits_files(dir: &Path) -> Vec<fs::DirEntry> {
    let mut entries: Vec<_> = fs::read_dir(dir)
        .expect("Failed to read directory")
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.path()
                .extension()
                .and_then(|ext| ext.to_str())
                .map(|ext| ext.eq_ignore_ascii_case("fits"))
                .unwrap_or(false)
        })
        .collect();
    entries.sort_by_key(|e| e.file_name());
    entries
}

fn load_images(
    entries: &[fs::DirEntry],
    crop: Option<f64>,
    config: &AppConfig,
    registered_stars: Option<&[RegisteredStar]>,
) -> Vec<ImageInfo> {
    let total = entries.len() as u64;
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{pos}/{len} [{bar:40.cyan/blue}] {msg}")
            .unwrap()
            .progress_chars("=>-"),
    );

    let mut images = Vec::new();
    let mut quality_sum = 0.0f64;
    let mut fwhm_sum = 0.0f64;
    let mut quality_count = 0usize;

    let mut last_images_window = vec![];
    let mut images_in_window = 100;
    let mut was_images_in_window_set = false;
    let mut last_quality = 0.0;

    for entry in entries {
        let file_path = entry.path();
        let file_name = entry.file_name().to_string_lossy().to_string();
        match AstroImage::load(&file_path, crop, config) {
            Ok(mut img) => {
                let independent_quality = img.quality();
                let independent_fwhm = img.fwhm();
                if !was_images_in_window_set {
                    images_in_window = (config.log_quality_window_t / img.exp_t()).ceil() as i32;
                    was_images_in_window_set = true
                }
                let constellation = registered_stars.map(|stars| {
                    Constellation::find_in_image_with_tolerance(
                        stars.to_vec(),
                        &img,
                        config.star_pattern_position_tolerance_px as f32,
                    )
                });
                let constellation_found = constellation.as_ref().map(|c| {
                    if c.found {
                        let quality_indices: HashSet<usize> = c
                            .registered_stars
                            .iter()
                            .enumerate()
                            .filter(|(_, rs)| rs.use_in_quality)
                            .filter_map(|(i, _)| c.star_mapping.get(&i).copied())
                            .collect();
                        img.recalculate_quality_for_star_indices(&quality_indices, config);
                        true
                    } else {
                        false
                    }
                });
                let matched_indices: Vec<usize> = constellation
                    .as_ref()
                    .filter(|c| c.found)
                    .map(|c| c.star_mapping.values().copied().collect())
                    .unwrap_or_default();
                let quality_indices: HashSet<usize> = constellation
                    .as_ref()
                    .filter(|c| c.found)
                    .map(|c| {
                        c.registered_stars
                            .iter()
                            .enumerate()
                            .filter(|(_, s)| s.use_in_quality)
                            .filter_map(|(i, _)| c.star_mapping.get(&i).copied())
                            .collect()
                    })
                    .unwrap_or_default();
                let matched_star_count = matched_indices.len();
                let star_brightness_adu = constellation
                    .as_ref()
                    .and_then(|c| star_pattern_brightness_adu(img.stars(), c));
                let background_raw_adu = img.background_raw_adu();
                let background_corrected_adu =
                    (background_raw_adu - config.background_bias_adu).max(0.0);
                let snr = star_pattern_snr(star_brightness_adu, img.background_noise_adu());
                if independent_quality.is_finite() {
                    quality_sum += independent_quality;
                    fwhm_sum += independent_fwhm;
                    quality_count += 1;
                    last_images_window.push(independent_quality);
                    if last_images_window.len() > images_in_window as usize {
                        last_images_window.remove(0);
                    }
                    last_quality =
                        last_images_window.iter().sum::<f64>() / last_images_window.len() as f64
                }
                let avg_q = quality_sum / quality_count as f64;
                let avg_f = fwhm_sum / quality_count as f64;
                let last_time = last_images_window.len() as f64 * img.exp_t() / 60.0;
                let msg = format!(
                    "avg FWHM: {:.2}  avg quality: {:.2}%, quality from last {:.1}min: {:.2}% ",
                    avg_f,
                    avg_q * 100.0,
                    last_time,
                    last_quality * 100.0
                );
                pb.set_message(msg);
                pb.inc(1);
                let stars = img.stars();
                let brightest_star_photons = stars.first().map_or(0.0, |s| s.magnitude);
                let star5_photons = stars.get(4).map_or(0.0, |s| s.magnitude);
                images.push(ImageInfo {
                    fwhm: independent_fwhm,
                    quality: independent_quality,
                    quality_image: constellation_found
                        .filter(|found| *found)
                        .and_then(|_| img.quality_for_star_indices(&quality_indices, config)),
                    star_count: img.star_count(),
                    file_name,
                    file_path,
                    constellation_found,
                    matched_star_count,
                    star_brightness_adu,
                    snr,
                    brightest_star_adu: stars.first().map_or(0.0, |s| s.magnitude_adu),
                    background_raw_adu,
                    background_corrected_adu,
                    brightest_star_photons,
                    star5_photons,
                });
            }
            Err(e) => {
                pb.println(format!("Error loading {}: {}", file_name, e));
                pb.inc(1);
            }
        }
    }
    pb.finish_with_message("Done");
    images
}

fn percentile(sorted: &[f64], p: usize) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() - 1) as f64 * p as f64 / 100.0).round() as usize;
    sorted[idx]
}

fn write_quality_map(
    dir: &Path,
    images: &[ImageInfo],
    low_star_threshold: Option<usize>,
    config: &AppConfig,
) {
    let map_path = dir.join("quality_map.txt");
    let mut map_file = fs::File::create(&map_path).expect("Failed to create quality_map.txt");

    let mut quality_vals: Vec<f64> = images.iter().map(|i| i.quality).collect();
    quality_vals.sort_by(|a, b| a.total_cmp(b));
    let mut quality_image_vals: Vec<f64> = images.iter().filter_map(|i| i.quality_image).collect();
    quality_image_vals.sort_by(|a, b| a.total_cmp(b));
    let mut fwhm_vals: Vec<f64> = images.iter().map(|i| i.fwhm).collect();
    fwhm_vals.sort_by(|a, b| b.total_cmp(a));

    let fn_width = images
        .iter()
        .map(|i| i.file_name.len())
        .max()
        .unwrap_or(8)
        .max(8);
    writeln!(
        map_file,
        "{:<fn_width$}  {:>9}  {:>9}  {:>9}  {:>6}  {:>12}  {:>10}  note",
        "filename",
        "fwhm",
        "quality",
        "qual_img",
        "stars",
        "brightest",
        "star5",
        fn_width = fn_width
    )
    .expect("Failed to write header");
    writeln!(
        map_file,
        "{:-<fn_width$}  {:>9}  {:->9}  {:->9}  {:->6}  {:->12}  {:->10}  {:->4}",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        fn_width = fn_width
    )
    .expect("Failed to write separator");

    let mut sorted: Vec<&ImageInfo> = images.iter().collect();
    sorted.sort_by(|a, b| b.quality.total_cmp(&a.quality));
    for img in &sorted {
        let note = if img.constellation_found == Some(false) {
            "no_constellation"
        } else {
            match low_star_threshold {
                Some(t) if img.star_count < t => "low_stars",
                _ => "",
            }
        };
        let quality_image_str = img
            .quality_image
            .map(|q| format!("{:>9.4}%", q * 100.0))
            .unwrap_or_else(|| format!("{:>9}", "-"));
        let brightest = if img.brightest_star_photons > 0.0 {
            format!("{:>12.0}", img.brightest_star_photons)
        } else {
            format!("{:>12}", "0")
        };
        let star5 = if img.star5_photons > 0.0 {
            format!("{:>10.0}", img.star5_photons)
        } else {
            format!("{:>10}", "0")
        };
        writeln!(
            map_file,
            "{:<fn_width$}  {:>9.4}  {:>9.4}%  {}  {:>6}  {}  {}  {}",
            img.file_name,
            img.fwhm,
            img.quality * 100.0,
            quality_image_str,
            img.star_count,
            brightest,
            star5,
            note,
            fn_width = fn_width
        )
        .expect("Failed to write quality map");
    }

    writeln!(map_file, "# Percentiles").expect("Failed to write");
    writeln!(
        map_file,
        "# {:>4}  {:>9}  {:>9}  {:>9}",
        "pct", "quality", "qual_img", "fwhm"
    )
    .expect("Failed to write");
    for p in (0..=100).step_by(5) {
        let q = percentile(&quality_vals, p);
        let qi = if quality_image_vals.is_empty() {
            String::from("         -")
        } else {
            format!("{:>9.4}%", percentile(&quality_image_vals, p) * 100.0)
        };
        let fwhm_p = percentile(&fwhm_vals, p);
        if p == 50 {
            println!("Median FWHM (SEEING): {:.4}", fwhm_from_quality(q));
            println!("MEDIAN QUALITY FOR SEQUENCE: {:.4} %", q * 100.0);
        }
        writeln!(
            map_file,
            "# {:>3}%  {:>9.4}%  {}  {:>9.4}",
            p,
            q * 100.0,
            qi,
            fwhm_p
        )
        .expect("Failed to write percentile");
    }
    writeln!(map_file, "#").expect("Failed to write");

    // Trend section: rolling averages sorted by filename (chronological proxy)
    let window = config.rolling_avg_window;
    if window > 0 && images.len() > window {
        let mut by_filename: Vec<&ImageInfo> = images.iter().collect();
        by_filename.sort_by(|a, b| a.file_name.cmp(&b.file_name));

        writeln!(map_file, "\n# Trend (rolling average, window = {})", window)
            .expect("Failed to write trend header");
        writeln!(map_file, "# images          quality   qual_img       fwhm")
            .expect("Failed to write trend header");

        let mut start = 0;
        while start < by_filename.len() {
            let end = (start + window).min(by_filename.len());
            let slice = &by_filename[start..end];
            let avg_q = slice.iter().map(|i| i.quality).sum::<f64>() / slice.len() as f64;
            let avg_fwhm = slice.iter().map(|i| i.fwhm).sum::<f64>() / slice.len() as f64;
            let qi_vals: Vec<f64> = slice.iter().filter_map(|i| i.quality_image).collect();
            let avg_qi_str = if qi_vals.is_empty() {
                format!("{:>9}", "-")
            } else {
                format!(
                    "{:>9.4}%",
                    qi_vals.iter().sum::<f64>() / qi_vals.len() as f64 * 100.0
                )
            };
            writeln!(
                map_file,
                "# {}-{}  {:>9.4}%  {}  {:>9.4}",
                start + 1,
                end,
                avg_q * 100.0,
                avg_qi_str,
                avg_fwhm
            )
            .expect("Failed to write trend row");
            start += window;
        }
    }

    writeln!(map_file, "#").expect("Failed to write");
    println!("\nQuality map written to: {}", map_path.display());
}

/* legacy selection helpers removed: filtering is implemented by metrics::FilterRule
and must not be confused with the removed --take/--remove modes. */

fn image_to_metrics(image: &ImageInfo) -> metrics::MetricValues {
    metrics::MetricValues {
        quality: image.quality,
        fwhm: image.fwhm,
        star_count: image.star_count,
        brightest_star_adu: image.brightest_star_adu,
        brightest_star_photons: image.brightest_star_photons,
        star5_photons: image.star5_photons,
        background_raw_adu: image.background_raw_adu,
        background_corrected_adu: image.background_corrected_adu,
        quality_star_pattern: image.quality_image,
        quality_star_pattern_source: image.quality_image.is_some(),
        star_brightness_adu: image.star_brightness_adu,
        snr: image.snr,
        star_pattern_found: image.constellation_found == Some(true),
        matched_star_count: image.matched_star_count,
    }
}

fn build_rules(args: &Args, has_pattern: bool) -> Result<Vec<metrics::FilterRule>, String> {
    use metrics::Metric::*;
    let mut rules = Vec::new();
    let mut add = |metric: metrics::Metric,
                   relative: Option<f64>,
                   absolute: Option<f64>|
     -> Result<(), String> {
        if relative.is_some() && absolute.is_some() {
            return Err(format!(
                "both relative and absolute thresholds supplied for {}",
                metric.key()
            ));
        }
        if let Some(v) = relative {
            rules.push(metrics::FilterRule::relative(metric, v)?);
        }
        if let Some(v) = absolute {
            rules.push(metrics::FilterRule::absolute(metric, v)?);
        }
        Ok(())
    };
    add(Snr, args.snr, args.snr_absolute)?;
    add(
        QualityStarPattern,
        args.quality_star_pattern,
        args.quality_star_pattern_absolute,
    )?;
    add(Background, args.background, args.background_absolute)?;
    add(
        StarBrightness,
        args.star_brightness,
        args.star_brightness_absolute,
    )?;
    add(Quality, args.quality, args.quality_absolute)?;
    add(Fwhm, args.fwhm, args.fwhm_absolute)?;
    if args.filter && rules.is_empty() {
        rules.push(metrics::FilterRule::relative(QualityStarPattern, 0.8)?);
        rules.push(metrics::FilterRule::relative(Background, 2.0)?);
        rules.push(metrics::FilterRule::relative(Snr, 0.7)?);
        rules.push(metrics::FilterRule::relative(StarBrightness, 0.7)?);
    }
    if rules.iter().any(|r| r.metric.requires_pattern()) && !has_pattern {
        return Err("filter requires a star pattern, but no pattern is available".into());
    }
    Ok(rules)
}

fn print_statistics(records: &[metrics::MetricRecord]) {
    println!("\n--- METRIC STATISTICS ({} records) ---", records.len());
    for metric in [
        metrics::Metric::Quality,
        metrics::Metric::Fwhm,
        metrics::Metric::QualityStarPattern,
        metrics::Metric::Background,
        metrics::Metric::StarBrightness,
        metrics::Metric::Snr,
    ] {
        let mut values: Vec<f64> = records
            .iter()
            .filter_map(|r| metric.value(&r.metrics))
            .filter(|v| v.is_finite())
            .collect();
        if values.is_empty() {
            continue;
        }
        let count = values.len();
        let mean = values.iter().sum::<f64>() / count as f64;
        values.sort_by(|a, b| a.total_cmp(b));
        println!(
            "{:<24} available={}/{} median={:.4} mean={:.4} min={:.4} max={:.4}",
            metric.key(),
            count,
            records.len(),
            metrics::median(&mut values.clone()).unwrap(),
            mean,
            values[0],
            values[count - 1]
        );
    }
    if records.iter().any(|r| !r.metrics.star_pattern_found) {
        println!(
            "Warning: images without a matched pattern have unavailable pattern-dependent metrics and will fail those filters."
        );
    }
}

fn bitmap_text(image: &mut image::RgbImage, x: u32, y: u32, text: &str, color: image::Rgb<u8>) {
    fn glyph(c: char) -> [u8; 7] {
        match c.to_ascii_uppercase() {
            'A' => [
                0b01110, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001,
            ],
            'B' => [
                0b11110, 0b10001, 0b11110, 0b10001, 0b10001, 0b10001, 0b11110,
            ],
            'C' => [
                0b01111, 0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b01111,
            ],
            'D' => [
                0b11110, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b11110,
            ],
            'E' => [
                0b11111, 0b10000, 0b11110, 0b10000, 0b10000, 0b10000, 0b11111,
            ],
            'F' => [
                0b11111, 0b10000, 0b11110, 0b10000, 0b10000, 0b10000, 0b10000,
            ],
            'G' => [
                0b01111, 0b10000, 0b10000, 0b10111, 0b10001, 0b10001, 0b01111,
            ],
            'H' => [
                0b10001, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001,
            ],
            'I' => [
                0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b11111,
            ],
            'K' => [
                0b10001, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010, 0b10001,
            ],
            'L' => [
                0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b11111,
            ],
            'M' => [
                0b10001, 0b11011, 0b10101, 0b10101, 0b10001, 0b10001, 0b10001,
            ],
            'N' => [
                0b10001, 0b11001, 0b10101, 0b10011, 0b10001, 0b10001, 0b10001,
            ],
            'O' => [
                0b01110, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110,
            ],
            'P' => [
                0b11110, 0b10001, 0b10001, 0b11110, 0b10000, 0b10000, 0b10000,
            ],
            'Q' => [
                0b01110, 0b10001, 0b10001, 0b10001, 0b10101, 0b10010, 0b01101,
            ],
            'R' => [
                0b11110, 0b10001, 0b11110, 0b10100, 0b10010, 0b10001, 0b10001,
            ],
            'S' => [
                0b01111, 0b10000, 0b10000, 0b01110, 0b00001, 0b00001, 0b11110,
            ],
            'T' => [
                0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100,
            ],
            'U' => [
                0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110,
            ],
            'W' => [
                0b10001, 0b10001, 0b10001, 0b10101, 0b10101, 0b11011, 0b10001,
            ],
            'Y' => [
                0b10001, 0b10001, 0b01010, 0b00100, 0b00100, 0b00100, 0b00100,
            ],
            '0' => [
                0b01110, 0b10011, 0b10101, 0b10101, 0b10101, 0b11001, 0b01110,
            ],
            '1' => [
                0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110,
            ],
            '2' => [
                0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111,
            ],
            '3' => [
                0b11110, 0b00001, 0b00001, 0b01110, 0b00001, 0b00001, 0b11110,
            ],
            '4' => [
                0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010,
            ],
            '5' => [
                0b11111, 0b10000, 0b10000, 0b11110, 0b00001, 0b00001, 0b11110,
            ],
            '6' => [
                0b00111, 0b01000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110,
            ],
            '7' => [
                0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b01000, 0b01000,
            ],
            '8' => [
                0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110,
            ],
            '9' => [
                0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00010, 0b11100,
            ],
            '.' => [0, 0, 0, 0, 0, 0, 0b00100],
            '_' => [0, 0, 0, 0, 0, 0, 0b11111],
            '-' => [0, 0, 0, 0b11111, 0, 0, 0],
            ' ' => [0; 7],
            _ => [0; 7],
        }
    }
    let mut cursor = x;
    for c in text.chars() {
        let rows = glyph(c);
        for (row, bits) in rows.iter().enumerate() {
            for col in 0..5 {
                if bits & (1 << (4 - col)) != 0
                    && cursor + col < image.width()
                    && y + (row as u32) < image.height()
                {
                    image.put_pixel(cursor + col, y + (row as u32), color);
                }
            }
        }
        cursor += 6;
    }
}

fn draw_metric_chart(
    dir: &Path,
    records: &[metrics::MetricRecord],
    metric: metrics::Metric,
) -> Result<(), String> {
    use image::{ImageBuffer, Rgb};
    let mut ordered: Vec<&metrics::MetricRecord> = records.iter().collect();
    ordered.sort_by(|a, b| {
        a.modified_ns
            .cmp(&b.modified_ns)
            .then_with(|| a.file_name.cmp(&b.file_name))
    });
    let values: Vec<Option<f64>> = ordered
        .iter()
        .map(|r| metric.value(&r.metrics).filter(|v| v.is_finite()))
        .collect();
    if values.iter().all(Option::is_none) {
        return Ok(());
    }
    let available: Vec<f64> = values.iter().filter_map(|v| *v).collect();
    let median = metrics::median(&mut available.clone()).unwrap_or(0.0);
    let suggested = match metric {
        metrics::Metric::QualityStarPattern => Some(median * 0.8),
        metrics::Metric::Background => Some(median * 2.0),
        metrics::Metric::StarBrightness => Some(median * 0.7),
        _ => None,
    };
    let mut scale_values = available.clone();
    if let Some(value) = suggested {
        scale_values.push(value);
    }
    let width = 900u32;
    let height = 420u32;
    let left = 55u32;
    let right = 20u32;
    let top = 35u32;
    let bottom = 35u32;
    let mut out: ImageBuffer<Rgb<u8>, Vec<u8>> =
        ImageBuffer::from_pixel(width, height, Rgb([255u8, 255u8, 255u8]));
    let min = scale_values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = scale_values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let range = (max - min).max(1e-12);
    for x in left..(width - right) {
        out.put_pixel(x, height - bottom, Rgb([0, 0, 0]));
    }
    for y in top..(height - bottom) {
        out.put_pixel(left, y, Rgb([0, 0, 0]));
    }
    // Five Y-axis ticks make the value scale readable even in the bitmap renderer.
    for tick in 0..=4 {
        let fraction = tick as f64 / 4.0;
        let y =
            height - bottom - 1 - (fraction * (height - top - bottom - 1) as f64).round() as u32;
        for dx in 0..5 {
            out.put_pixel(left.saturating_sub(dx), y, Rgb([0, 0, 0]));
        }
        let value = min + fraction * range;
        bitmap_text(
            &mut out,
            2,
            y.saturating_sub(3),
            &format!("{value:.2}"),
            Rgb([0, 0, 0]),
        );
    }
    for (i, value) in values.iter().enumerate() {
        let x = left + (i as u32 * (width - left - right - 1) / values.len().max(1) as u32);
        let Some(value) = value else {
            continue;
        }; // missing values are gaps
        let y = height
            - bottom
            - 1
            - (((*value - min) / range) * (height - top - bottom - 1) as f64)
                .round()
                .clamp(0.0, (height - top - bottom - 1) as f64) as u32;
        for dx in 0..5 {
            for dy in 0..5 {
                if x + dx < width && y + dy < height {
                    out.put_pixel(x + dx, y + dy, Rgb([20, 70, 180]));
                }
            }
        }
    }
    bitmap_text(&mut out, left + 5, 5, metric.key(), Rgb([0, 0, 0]));
    bitmap_text(
        &mut out,
        width - 300,
        5,
        &format!("M MEDIAN {median:.2}"),
        Rgb([220, 0, 0]),
    );
    if let Some(value) = suggested {
        bitmap_text(
            &mut out,
            width - 300,
            15,
            &format!("T THRESHOLD {value:.2}"),
            Rgb([0, 150, 0]),
        );
    }
    let line = |out: &mut ImageBuffer<Rgb<u8>, Vec<u8>>, value: f64, color: Rgb<u8>| {
        let y = height
            - bottom
            - 1
            - (((value - min) / range) * (height - top - bottom - 1) as f64)
                .round()
                .clamp(0.0, (height - top - bottom - 1) as f64) as u32;
        for x in left..(width - right) {
            out.put_pixel(x, y, color);
        }
    };
    line(&mut out, median, Rgb([220, 0, 0]));
    if let Some(value) = suggested {
        line(&mut out, value, Rgb([0, 150, 0]));
    }
    out.save(dir.join(format!("metrics_{}.png", metric.key())))
        .map_err(|e| e.to_string())
}

fn filter_files(
    dir: &Path,
    records: &[metrics::MetricRecord],
    rules: &[metrics::FilterRule],
) -> Result<(), String> {
    let medians = metrics::medians(records);
    for rule in rules {
        let threshold = rule.threshold(medians.get(&rule.metric).copied())?;
        let median = medians.get(&rule.metric).copied().unwrap_or(0.0);
        let direction = if rule.metric.higher_is_better() {
            ">="
        } else {
            "<="
        };
        if let Some(multiplier) = rule.relative {
            println!(
                "{}: median={:.4}, multiplier={:.3}, threshold={:.4}, direction keep {} threshold",
                rule.metric.key(),
                median,
                multiplier,
                threshold,
                direction
            );
        } else {
            println!(
                "{}: absolute={:.4}, threshold={:.4}, direction keep {} threshold",
                rule.metric.key(),
                rule.absolute.unwrap_or(threshold),
                threshold,
                direction
            );
        }
    }
    let rejected: Vec<_> = records
        .iter()
        .filter(|r| !metrics::passes_all_filters(r, rules, &medians).unwrap_or(false))
        .collect();
    let destination = metrics::unique_removed_folder(dir, rules, metrics::current_timestamp());
    if rejected.is_empty() {
        println!(
            "Kept {} files; no files were rejected or moved.",
            records.len()
        );
        return Ok(());
    }
    fs::create_dir_all(&destination).map_err(|e| e.to_string())?;
    for record in &rejected {
        let source = dir.join(&record.file_name);
        let target = destination.join(&record.file_name);
        fs::rename(&source, &target)
            .or_else(|_| {
                fs::copy(&source, &target)
                    .map(|_| ())
                    .and_then(|_| fs::remove_file(&source))
            })
            .map_err(|e| e.to_string())?;
    }
    println!(
        "Kept {}, moved {} files to {}",
        records.len() - rejected.len(),
        rejected.len(),
        destination.display()
    );
    Ok(())
}

/// Parse the two arguments of --divide.  Keeping this separate from clap also
/// makes the accepted metric names and validation explicit to the user.
fn parse_divide(value: &[String]) -> Result<(metrics::Metric, usize), String> {
    if value.len() != 2 {
        return Err("--divide requires METRIC and COUNT, e.g. --divide snr 1000".into());
    }
    let metric_name = value[0].to_ascii_lowercase().replace('-', "_");
    let metric = match metric_name.as_str() {
        "quality" => metrics::Metric::Quality,
        "fwhm" => metrics::Metric::Fwhm,
        "quality_star_pattern" => metrics::Metric::QualityStarPattern,
        "background" => metrics::Metric::Background,
        "star_brightness" => metrics::Metric::StarBrightness,
        "snr" => metrics::Metric::Snr,
        _ => return Err(format!(
            "unknown --divide metric '{}'; choose quality, fwhm, quality_star_pattern, background, star_brightness or snr",
            value[0]
        )),
    };
    let count = value[1].parse::<usize>().map_err(|_| {
        format!("invalid --divide COUNT '{}'; COUNT must be a positive integer", value[1])
    })?;
    if count == 0 {
        return Err("--divide COUNT must be greater than zero".into());
    }
    Ok((metric, count))
}

fn divide_files(
    dir: &Path,
    records: &[metrics::MetricRecord],
    metric: metrics::Metric,
    batch_size: usize,
) -> Result<(), String> {
    // At this point rejected files have already been moved.  Checking the
    // source path makes this function operate only on files which passed.
    let mut kept: Vec<&metrics::MetricRecord> = records
        .iter()
        .filter(|record| dir.join(&record.file_name).is_file())
        .collect();
    kept.sort_by(|a, b| {
        let av = metric.value(&a.metrics).filter(|v| v.is_finite());
        let bv = metric.value(&b.metrics).filter(|v| v.is_finite());
        match (av, bv) {
            (Some(a_value), Some(b_value)) => {
                let ordering = if metric.higher_is_better() { b_value.total_cmp(&a_value) } else { a_value.total_cmp(&b_value) };
                ordering.then_with(|| a.file_name.cmp(&b.file_name))
            }
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => a.file_name.cmp(&b.file_name),
        }
    });

    if kept.iter().any(|r| metric.value(&r.metrics).filter(|v| v.is_finite()).is_none()) {
        return Err(format!(
            "cannot divide files: metric '{}' is unavailable for one or more kept files",
            metric.key()
        ));
    }
    for (batch_index, chunk) in kept.chunks(batch_size).enumerate() {
        let folder = dir.join(format!("{}_{}", batch_index + 1, metric.key()));
        fs::create_dir_all(&folder).map_err(|e| e.to_string())?;
        for record in chunk {
            let source = dir.join(&record.file_name);
            let target = folder.join(&record.file_name);
            fs::rename(&source, &target)
                .or_else(|_| fs::copy(&source, &target).map(|_| ()).and_then(|_| fs::remove_file(&source)))
                .map_err(|e| e.to_string())?;
        }
        println!("{}: moved {} files to {}", folder.file_name().unwrap().to_string_lossy(), chunk.len(), folder.display());
    }
    Ok(())
}

fn process_directory(
    dir: &Path,
    args: &Args,
    registered_stars: Option<&[RegisteredStar]>,
    config: &AppConfig,
) -> Result<(), String> {
    if args.divide.is_some() && !args.filter {
        return Err("--divide must be used together with --filter".into());
    }
    let divide = args.divide.as_deref().map(parse_divide).transpose()?;
    let filter_rules = if args.filter {
        match build_rules(args, registered_stars.is_some()) {
            Ok(rules) => Some(rules),
            Err(e) => {
                return Err(e);
            }
        }
    } else {
        None
    };
    let all_entries = collect_fits_files(dir);
    println!("Input files discovered: {}", all_entries.len());
    if all_entries.is_empty() {
        return Err("no FITS files found".into());
    }
    let entries: Vec<_> = if let Some(n) = args.check_count {
        all_entries.into_iter().rev().take(n).collect()
    } else { all_entries };
    println!("Analysis started");
    let records: Vec<metrics::MetricRecord> = load_images(&entries, args.crop, config, registered_stars)
        .iter().map(|image| metrics::MetricRecord {
            file_name: image.file_name.clone(),
            modified_ns: fs::metadata(&image.file_path).ok().and_then(|m| m.modified().ok()).and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok()).map(|d| d.as_nanos()).unwrap_or(0),
            metrics: image_to_metrics(image),
        }).collect();
    let images: Vec<_> = records.iter().map(|r| ImageInfo {
        file_name: r.file_name.clone(), file_path: dir.join(&r.file_name), quality: r.metrics.quality, fwhm: r.metrics.fwhm,
        quality_image: r.metrics.quality_star_pattern.filter(|_| r.metrics.quality_star_pattern_source), star_count: r.metrics.star_count,
        constellation_found: Some(r.metrics.star_pattern_found), matched_star_count: r.metrics.matched_star_count, star_brightness_adu: r.metrics.star_brightness_adu, snr: r.metrics.snr,
        brightest_star_adu: r.metrics.brightest_star_adu, background_raw_adu: r.metrics.background_raw_adu, background_corrected_adu: r.metrics.background_corrected_adu,
        brightest_star_photons: r.metrics.brightest_star_photons, star5_photons: r.metrics.star5_photons }).collect();
    if args.check_count.is_none() { write_quality_map(dir, &images, None, config); }
    print_statistics(&records);
    println!("Analysis finished");
    if args.check_count.is_none()
        && !args.filter
        && build_rules(args, registered_stars.is_some())
            .map(|r| !r.is_empty())
            .unwrap_or(false)
    {
        eprintln!("Warning: thresholds supplied without --filter; no files moved.");
    }
    if args.check_count.is_none()
        && args.filter
        && let Some(rules) = filter_rules
    {
        filter_files(dir, &records, &rules)?;
        if let Some((metric, batch_size)) = divide {
            divide_files(dir, &records, metric, batch_size)?;
        }
    }
    if args.check_count.is_none() {
        // Filtering may have moved files, so chart only the files still in the directory.
        let chart_records: Vec<_> = records
            .iter()
            .filter(|r| dir.join(&r.file_name).is_file())
            .cloned()
            .collect();
        for metric in [
            metrics::Metric::Quality,
            metrics::Metric::Fwhm,
            metrics::Metric::QualityStarPattern,
            metrics::Metric::Background,
            metrics::Metric::StarBrightness,
            metrics::Metric::Snr,
        ] {
            let _ = draw_metric_chart(dir, &chart_records, metric);
        }
    }
    Ok(())
}

fn make_star_pattern(
    folder: &Path,
    config: &AppConfig,
    explicit_numbers: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    let entries = collect_fits_files(folder);
    if entries.is_empty() {
        return Err("folder contains no FITS files".into());
    }
    let sample_paths =
        star_pattern::select_sample_paths(entries.iter().map(|e| e.path()).collect(), 10, 0);
    let mut loaded = Vec::new();
    for path in sample_paths {
        match AstroImage::load(&path, None, config) {
            Ok(img) => loaded.push(img),
            Err(e) => eprintln!("Warning: skipping {}: {}", path.display(), e),
        }
    }
    let reference = loaded.first().ok_or("unable to load FITS samples")?;
    let to_sample = |img: &AstroImage| {
        img.stars()
            .iter()
            .map(|s| star_pattern::PatternStarSample {
                x: s.pos.x,
                y: s.pos.y,
                magnitude: s.magnitude,
                magnitude_adu: s.magnitude_adu,
                // Detector-scale brightness, as required by the 10..50% preference.
                brightest_pixel_part: (s.brightest_pixel_adu / u16::MAX as f64).clamp(0.0, 1.0),
            })
            .collect::<Vec<_>>()
    };
    let frames: Vec<Vec<_>> = loaded.iter().map(to_sample).collect();
    let reference_stars = frames.first().cloned().unwrap_or_default();
    let aggregated = star_pattern::aggregate_samples(
        &reference_stars,
        &frames,
        config.star_pattern_position_tolerance_px,
    );
    let candidates: Vec<_> = aggregated.iter().map(|a| a.sample.clone()).collect();
    let rec = star_pattern::recommend_stars(&candidates, reference.width(), reference.height());
    if aggregated.len() < 2 {
        let _ = reference.save_stars_jpg(folder.join("star_pattern_reference_problem.jpg"));
        return Err(format!("unable to generate pattern: fewer than 2 stars remain in all sampled frames ({}).", aggregated.len()).into());
    }
    if rec.len() < 3 {
        println!("Only {} stable stars available; recommendation requires 3.", aggregated.len());
    }
    let recommended_numbers: Vec<usize> = rec
        .iter()
        .map(|&i| aggregated[i].reference_index + 1)
        .collect();
    reference.save_stars_jpg(folder.join("star_pattern_candidates.jpg"))?;
    println!(
        "Recommended stars (numbers in star_pattern_candidates.jpg): {:?}",
        recommended_numbers
    );

    let numbers = if let Some(raw) = explicit_numbers {
        raw.to_string()
    } else {
        println!("Accept recommendation? [Y/n] (or enter comma-separated JPEG star numbers)");
        let mut answer = String::new();
        io::stdin().read_line(&mut answer)?;
        if answer.trim().is_empty() || answer.trim().eq_ignore_ascii_case("y") {
            recommended_numbers
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        } else if answer.trim().eq_ignore_ascii_case("n") {
            println!(
                "Enter at least two star numbers from star_pattern_candidates.jpg (for example 1,4):"
            );
            let mut manual = String::new();
            io::stdin().read_line(&mut manual)?;
            manual
        } else {
            answer
        }
    };
    let selected: Vec<usize> = numbers
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<usize>())
        .collect::<Result<_, _>>()?;
    if selected.len() < 2 {
        return Err("a valid star pattern requires at least 2 selected stars".into());
    }
    if selected
            .iter()
            .any(|n| *n == 0 || *n > reference_stars.len())
    {
        return Err(format!(
            "star number must be between 1 and {} in the candidate JPEG",
            reference_stars.len()
        )
            .into());
    }
    let mut pattern = Vec::new();
    for number in selected {
        let reference_index = number - 1;
        let aggregate = aggregated
            .iter()
            .find(|a| a.reference_index == reference_index)
            .ok_or_else(|| {
                format!(
                    "star {} is a valid JPEG number but is unstable across sampled frames",
                    number
                )
            })?;
        let s = &aggregate.sample;
        pattern.push(star_pattern::StarPatternEntry {
            x: s.x,
            y: s.y,
            magnitude: s.magnitude,
            use_in_quality: true,
            median_brightness_adu: Some(s.magnitude_adu),
            median_brightest_pixel_part: Some(s.brightest_pixel_part),
        });
    }
    star_pattern::validate_pattern(&pattern, config.star_pattern_position_tolerance_px)?;
    let output = folder.join("stars_pattern.json");
    fs::write(&output, serde_json::to_string_pretty(&pattern)?)?;
    println!("Pattern written to {}", output.display());
    Ok(())
}

fn main() {
    println!("Starting lucky-star");
    let args = Args::parse();
    println!("Arguments parsed");

    // Pattern generation is independent of any pre-existing pattern file.
    if let Some(folder) = &args.make_star_pattern {
        let config = AppConfig::load_or_default();
    println!("Configuration loaded");
        if let Err(e) = config.validate() {
            eprintln!("Error: {}", e);
            std::process::exit(2);
        }
        if let Err(e) = make_star_pattern(
            Path::new(folder),
            &config,
            args.star_pattern_numbers.as_deref(),
        ) {
            eprintln!("Error creating star pattern: {}", e);
            std::process::exit(1);
        }
        return;
    }

    let pattern_file: Option<PathBuf> = args
        .star_pattern
        .as_ref()
        .map(PathBuf::from)
        .or_else(|| star_pattern::default_pattern_path(Path::new(&args.path)));
    let registered_stars: Option<Vec<RegisteredStar>> = match &pattern_file {
        Some(p) => Some(load_stars_from_json(p).unwrap_or_else(|e| {
            eprintln!("Error loading star pattern '{}': {}", p.display(), e);
            std::process::exit(1);
        })),
        None => {
            let default_stars = star_pattern::default_pattern_path(Path::new(&args.path));
            if let Some(default_stars) = default_stars {
                match load_stars_from_json(&default_stars) {
                    Ok(stars) => {
                        println!("Using star pattern {}.", default_stars.display());
                        Some(stars)
                    }
                    Err(e) => {
                        eprintln!("Warning: failed to load stars.json: {}", e);
                        None
                    }
                }
            } else {
                None
            }
        }
    };

    let config = AppConfig::load_or_default();
    if let Err(e) = config.validate() {
        eprintln!("Error: {}", e);
        std::process::exit(2);
    }
    let path = Path::new(&args.path);
    if path.is_file() {
        process_single_file(
            path,
            args.crop,
            args.save_stars,
            &config,
            registered_stars.as_deref(),
        );
    } else if path.is_dir() {
        if args.save_stars {
            let entries = collect_fits_files(path);
            for entry in entries {
                process_single_file(
                    &entry.path(),
                    args.crop,
                    true,
                    &config,
                    registered_stars.as_deref(),
                );
            }
        } else {
            if let Err(e) = process_directory(path, &args, registered_stars.as_deref(), &config) {
                eprintln!("Error: {}", e);
                std::process::exit(1);
            }
        }
    } else {
        eprintln!("Error: path does not exist: {}", args.path);
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use vector2d::Vector2D;

    #[test]
    fn filter_does_not_create_removed_folder_when_every_record_passes() {
        let dir = std::env::temp_dir().join(format!(
            "lucky-star-filter-test-{}-{}",
            std::process::id(),
            metrics::current_timestamp()
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("a.fits"), b"fixture").unwrap();
        let record = metrics::MetricRecord {
            file_name: "a.fits".into(),
            modified_ns: 1,
            metrics: metrics::MetricValues {
                quality: 1.0,
                fwhm: 1.0,
                star_count: 1,
                brightest_star_adu: 1.0,
                brightest_star_photons: 1.0,
                star5_photons: 0.0,
                background_raw_adu: 1.0,
                background_corrected_adu: 1.0,
                quality_star_pattern: None,
                quality_star_pattern_source: false,
                star_brightness_adu: None,
                snr: None,
                star_pattern_found: false,
                matched_star_count: 0,
            },
        };
        let rule = metrics::FilterRule::absolute(metrics::Metric::Quality, 0.5).unwrap();
        filter_files(&dir, &[record], &[rule]).unwrap();
        assert!(!dir.join("removed_quality_0.500").exists());
        assert!(dir.join("a.fits").exists());
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn disabled_pattern_stars_do_not_affect_brightness_or_snr() {
        let stars = vec![
            Star {
                pos: Vector2D::new(1, 1),
                magnitude: 1.0,
                magnitude_adu: 100.0,
                brightest_pixel_adu: 1.0,
                brightest_pixel_part: 1.0,
                top_4_pixels_part: 1.0,
                ill_defined: false,
            },
            Star {
                pos: Vector2D::new(2, 2),
                magnitude: 1.0,
                magnitude_adu: 900.0,
                brightest_pixel_adu: 1.0,
                brightest_pixel_part: 1.0,
                top_4_pixels_part: 1.0,
                ill_defined: false,
            },
        ];
        let constellation = Constellation {
            registered_stars: vec![
                RegisteredStar {
                    pos: Vector2D::new(1, 1),
                    magnitude: 1.0,
                    use_in_quality: true,
                    median_brightness_adu: Some(100.0),
                    median_brightest_pixel_part: None,
                },
                RegisteredStar {
                    pos: Vector2D::new(2, 2),
                    magnitude: 1.0,
                    use_in_quality: false,
                    median_brightness_adu: Some(900.0),
                    median_brightest_pixel_part: None,
                },
            ],
            found: true,
            star_mapping: HashMap::from([(0, 0), (1, 1)]),
            transform: None,
            position_tolerance_px: 7.0,
        };
        let signal = star_pattern_brightness_adu(&stars, &constellation);
        assert_eq!(signal, Some(100.0));
        assert_eq!(star_pattern_snr(signal, 25.0), Some(4.0));
    }
}

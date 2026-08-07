use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use serde::{Deserialize, Serialize};
use crate::astro_image::{fwhm_from_quality, AstroImage};
use crate::constellation::{Constellation, RegisteredStar, load_stars_from_json};
use crate::helpers::median;

pub mod astro_image;
pub mod star;
pub mod helpers;
pub mod constellation;
pub mod star_pattern;
pub mod metrics;


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
            background_bias_adu: 0.0,
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
            star_pattern_position_tolerance_px: partial.star_pattern_position_tolerance_px.unwrap_or(default.star_pattern_position_tolerance_px),
            background_bias_adu: partial.background_bias_adu.unwrap_or(default.background_bias_adu),
        })
    }
}

impl AppConfig {
    fn validate(&self) -> Result<(), String> {
        if self.background_bias_adu < 0.0 { return Err("background_bias_adu must not be negative".into()); }
        if self.star_pattern_position_tolerance_px <= 0.0 { return Err("star_pattern_position_tolerance_px must be greater than zero".into()); }
        Ok(())
    }

    fn load_or_default() -> Self {
        // Try executable directory first
        if let Ok(exe_path) = std::env::current_exe() {
            if let Some(exe_dir) = exe_path.parent() {
                let config_path = exe_dir.join("config.json");
                if let Ok(content) = fs::read_to_string(&config_path) {
                    match serde_json::from_str(&content) {
                        Ok(cfg) => return cfg,
                        Err(e) => eprintln!("Warning: failed to parse config.json in executable directory: {}", e),
                    }
                }
            }
        }
        // Try current working directory
        if let Ok(cwd) = std::env::current_dir() {
            let config_path = cwd.join("config.json");
            if let Ok(content) = fs::read_to_string(&config_path) {
                match serde_json::from_str(&content) {
                    Ok(cfg) => return cfg,
                    Err(e) => eprintln!("Warning: failed to parse config.json in current directory: {}", e),
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

    #[arg(long)] snr: Option<f64>,
    #[arg(long)] quality_star_pattern: Option<f64>,
    #[arg(long)] background: Option<f64>,
    #[arg(long)] star_brightness: Option<f64>,
    #[arg(long)] quality: Option<f64>,
    #[arg(long)] fwhm: Option<f64>,
    #[arg(long)] snr_absolute: Option<f64>,
    #[arg(long)] background_absolute: Option<f64>,
    #[arg(long)] star_brightness_absolute: Option<f64>,
    #[arg(long)] quality_star_pattern_absolute: Option<f64>,
    #[arg(long)] quality_absolute: Option<f64>,
    #[arg(long)] fwhm_absolute: Option<f64>,

    /// Only search for stars in the central fraction of the image (e.g. 0.3 = central 30% width and height)
    #[arg(long, value_name = "FRACTION")]
    crop: Option<f64>,

    /// Save annotated star image to a JPG file
    #[arg(long, short)]
    save_stars: bool,

    /// Path to a JSON file with reference stars for constellation-based quality filtering
    #[arg(long, value_name = "FILE")]
    star_pattern: Option<String>,

    /// Quick analysis of the N newest FITS images; updates cache and skips charts.
    #[arg(short = 'c', long = "check-count", value_name = "N")]
    check_count: Option<usize>,

    /// Preserve the legacy seeing summary mode.
    #[arg(long, value_name = "N")]
    check_seeing: Option<usize>,

    /// Interactively generate stars_pattern.json from a session folder
    #[arg(long, value_name = "FOLDER")]
    make_star_pattern: Option<String>,
}

struct ImageInfo {
    file_name: String,
    file_path: PathBuf,
    modified_ns: u128,
    size: u64,
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



fn apply_constellation_quality(img: &mut AstroImage, registered_stars: &[RegisteredStar], config: &AppConfig, label: &str) -> bool {
    let constellation = Constellation::find_in_image_with_tolerance(registered_stars.to_vec(), img, config.star_pattern_position_tolerance_px as f32);
    if constellation.found {
        let quality_indices: HashSet<usize> = constellation.registered_stars.iter()
            .enumerate()
            .filter(|(_, rs)| rs.use_in_quality)
            .filter_map(|(i, _)| constellation.star_mapping.get(&i).copied())
            .collect();
        img.recalculate_quality_for_star_indices(&quality_indices, config);
        true
    } else {
        eprintln!("Warning: constellation not found in '{}', falling back to regular quality.", label);
        false
    }
}

fn process_single_file(path: &Path, crop: Option<f64>, save_stars: bool, config: &AppConfig, registered_stars: Option<&Vec<RegisteredStar>>) {
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
            e.path().extension()
                .and_then(|ext| ext.to_str())
                .map(|ext| ext.eq_ignore_ascii_case("fits"))
                .unwrap_or(false)
        })
        .collect();
    entries.sort_by_key(|e| e.file_name());
    entries
}

fn load_images(entries: &[fs::DirEntry], crop: Option<f64>, config: &AppConfig, registered_stars: Option<&Vec<RegisteredStar>>) -> Vec<ImageInfo> {
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
                let metadata = entry.metadata().ok();
                if !was_images_in_window_set {
                    images_in_window = (config.log_quality_window_t / img.exp_t()).ceil() as i32;
                    was_images_in_window_set = true
                }
                let constellation = registered_stars.map(|stars| Constellation::find_in_image_with_tolerance(stars.to_vec(), &img, config.star_pattern_position_tolerance_px as f32));
                let constellation_found = constellation.as_ref().map(|c| {
                    if c.found {
                        let quality_indices: HashSet<usize> = c.registered_stars.iter().enumerate()
                            .filter(|(_, rs)| rs.use_in_quality)
                            .filter_map(|(i, _)| c.star_mapping.get(&i).copied())
                            .collect();
                        img.recalculate_quality_for_star_indices(&quality_indices, config);
                        true
                    } else {
                        false
                    }
                });
                let matched_indices: Vec<usize> = constellation.as_ref().filter(|c| c.found).map(|c| c.star_mapping.values().copied().collect()).unwrap_or_default();
                let quality_indices: HashSet<usize> = constellation.as_ref().filter(|c| c.found).map(|c| c.registered_stars.iter().enumerate().filter(|(_, s)| s.use_in_quality).filter_map(|(i, _)| c.star_mapping.get(&i).copied()).collect()).unwrap_or_default();
                let matched_star_count = matched_indices.len();
                let star_brightness_adu = constellation_found.filter(|found| *found).map(|_| matched_indices.iter().filter_map(|i| img.stars().get(*i)).map(|s| s.magnitude_adu).sum::<f64>());
                let background_raw_adu = img.background_raw_adu();
                let background_corrected_adu = (background_raw_adu - config.background_bias_adu).max(0.0);
                let snr = star_brightness_adu.and_then(|signal| (background_corrected_adu > 0.0).then_some(signal / background_corrected_adu.sqrt()));
                if independent_quality.is_finite(){
                    quality_sum += independent_quality;
                    fwhm_sum += independent_fwhm;
                    quality_count += 1;
                    last_images_window.push(independent_quality);
                    if last_images_window.len() > images_in_window as usize {
                        last_images_window.remove(0);
                    }
                    last_quality = last_images_window.iter().sum::<f64>() / last_images_window.len() as f64
                }
                let avg_q = quality_sum / quality_count as f64;
                let avg_f = fwhm_sum / quality_count as f64;
                let last_time = last_images_window.len() as f64 * img.exp_t() / 60.0;
                let msg = format!("avg FWHM: {:.2}  avg quality: {:.2}%, quality from last {:.1}min: {:.2}% ",avg_f, avg_q * 100.0, last_time, last_quality * 100.0);
                pb.set_message(msg);
                pb.inc(1);
                let stars = img.stars();
                let brightest_star_photons = stars.first().map_or(0.0, |s| s.magnitude);
                let star5_photons = stars.get(4).map_or(0.0, |s| s.magnitude);
                images.push(ImageInfo {
                    fwhm: independent_fwhm,
                    quality: independent_quality,
                    quality_image: constellation_found.filter(|found| *found).and_then(|_| img.quality_for_star_indices(&quality_indices, config)),
                    star_count: img.star_count(),
                    file_name,
                    file_path,
                    size: metadata.as_ref().map(|m| m.len()).unwrap_or(0),
                    modified_ns: metadata.as_ref().map(metrics::modified_ns).unwrap_or(0),
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

fn write_quality_map(dir: &Path, images: &[ImageInfo], low_star_threshold: Option<usize>, config: &AppConfig) {
    let map_path = dir.join("quality_map.txt");
    let mut map_file = fs::File::create(&map_path).expect("Failed to create quality_map.txt");

    let mut quality_vals: Vec<f64> = images.iter().map(|i| i.quality).collect();
    quality_vals.sort_by(|a, b| a.total_cmp(b));
    let mut quality_image_vals: Vec<f64> = images.iter()
        .filter_map(|i| i.quality_image)
        .collect();
    quality_image_vals.sort_by(|a, b| a.total_cmp(b));
    let mut fwhm_vals: Vec<f64> = images.iter().map(|i| i.fwhm).collect();
    fwhm_vals.sort_by(|a, b| b.total_cmp(a));

    let fn_width = images.iter().map(|i| i.file_name.len()).max().unwrap_or(8).max(8);
    writeln!(
        map_file,
        "{:<fn_width$}  {:>9}  {:>9}  {:>9}  {:>6}  {:>12}  {:>10}  {}",
        "filename","fwhm", "quality", "qual_img", "stars", "brightest", "star5", "note",
        fn_width = fn_width
    ).expect("Failed to write header");
    writeln!(
        map_file,
        "{:-<fn_width$}  {:>9}  {:->9}  {:->9}  {:->6}  {:->12}  {:->10}  {:->4}",
        "","", "", "", "", "", "", "",
        fn_width = fn_width
    ).expect("Failed to write separator");

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
        let quality_image_str = img.quality_image
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
            img.file_name, img.fwhm, img.quality * 100.0, quality_image_str,
            img.star_count, brightest, star5, note,
            fn_width = fn_width
        )
        .expect("Failed to write quality map");
    }

    writeln!(map_file, "# Percentiles").expect("Failed to write");
    writeln!(map_file, "# {:>4}  {:>9}  {:>9}  {:>9}", "pct", "quality", "qual_img", "fwhm").expect("Failed to write");
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
        writeln!(map_file, "# {:>3}%  {:>9.4}%  {}  {:>9.4}", p, q * 100.0, qi, fwhm_p).expect("Failed to write percentile");
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
                format!("{:>9.4}%", qi_vals.iter().sum::<f64>() / qi_vals.len() as f64 * 100.0)
            };
            writeln!(map_file, "# {}-{}  {:>9.4}%  {}  {:>9.4}", start + 1, end, avg_q * 100.0, avg_qi_str, avg_fwhm)
                .expect("Failed to write trend row");
            start += window;
        }
    }

    writeln!(map_file, "#").expect("Failed to write");
    println!("\nQuality map written to: {}", map_path.display());
}

fn compute_star_threshold(images: &[ImageInfo]) -> (usize, usize) {
    let mut sorted_stars: Vec<usize> = images.iter().map(|i| i.star_count).collect();
    sorted_stars.sort_unstable();
    let median = median(&sorted_stars).unwrap_or_default();
    let threshold = (median as f64 * 0.7) as usize;
    (median, threshold)
}
fn select_best_by_percent(images: &[ImageInfo], take_pct: f64, median_stars: usize, low_star_threshold: usize, use_constellation: bool, keep_low: bool) -> HashSet<&str> {
    let take_pct = take_pct.clamp(0.0, 1.0);
    let total = images.len();
    let count_to_take = ((total as f64 * take_pct).ceil() as usize)
        .max(1)
        .min(total);

    let mut eligible: Vec<&ImageInfo> = if keep_low {
        images.iter().map(|i| i).collect()
    } else if use_constellation {
        images.iter().filter(|i| i.constellation_found == Some(true)).collect()
    } else {
        images.iter().filter(|i| i.star_count >= low_star_threshold).collect()
    };

    eligible.sort_by(|a, b| b.quality.total_cmp(&a.quality));

    let selected: HashSet<&str> = eligible.iter()
        .take(count_to_take)
        .map(|i| i.file_name.as_str())
        .collect();

    if use_constellation {
        println!(
            "Selection: {}/{} with constellation found, copying top {} by quality.",
            eligible.len(), total, selected.len()
        );
    } else {
        println!(
            "Selection: {}/{} eligible (median {} stars, min threshold {}), copying top {} by quality.",
            eligible.len(), total, median_stars, low_star_threshold, selected.len()
        );
    }

    selected
}

fn select_best_by_quality(images: &[ImageInfo], take_quality: f64, use_constellation: bool) -> HashSet<&str> {
    let take_quality = take_quality.clamp(0.0, 1.0);

    let mut eligible: Vec<&ImageInfo> = if use_constellation {
        images.iter().filter(|i| i.constellation_found == Some(true)).collect()
    } else {
        images.iter().collect()
    };

    eligible.sort_by(|a, b| {
        b.quality_image.unwrap_or(b.quality)
            .total_cmp(&a.quality_image.unwrap_or(a.quality))
    });

    let selected: HashSet<&str> = eligible.iter()
        .filter(|e|e.quality > take_quality)
        .map(|i| i.file_name.as_str())
        .collect();

    println!(
        "Quality-image selection: {}/{} eligible, copying top {} by quality_image.",
        eligible.len(), images.len(), selected.len()
    );

    selected
}

fn copy_to_named_folder(dir: &Path, images: &[ImageInfo], selected: &HashSet<&str>, folder_name: &str) {
    let dest_dir = dir.join(folder_name);
    fs::create_dir_all(&dest_dir).expect("Failed to create destination directory");
    let mut copied = 0usize;
    for img in images {
        if selected.contains(img.file_name.as_str()) {
            let dest = dest_dir.join(&img.file_name);
            if dest.exists() {
                continue
            }
            fs::copy(&img.file_path, &dest).expect("Failed to copy file");
            copied += 1;
        }
    }
    println!("Copied {} images to: {}", copied, dest_dir.display());
}

fn remove_original_images(dir: &Path, images: &[ImageInfo], selected: &HashSet<&str>) {
    let remove_dir = dir.join("remove");
    fs::create_dir_all(&remove_dir).expect("Failed to create 'remove' directory");
    let mut moved = 0usize;
    for img in images {
        if selected.contains(img.file_name.as_str()) {
            continue; // Keep selected originals in place
        }
        let dest = remove_dir.join(&img.file_name);
        if let Err(_) = fs::rename(&img.file_path, &dest) {
            // Cross-device fallback: copy then delete
            fs::copy(&img.file_path, &dest)
                .expect("Failed to copy file to remove directory");
            fs::remove_file(&img.file_path)
                .expect("Failed to remove original after copy");
        }
        moved += 1;
    }
    println!("Moved {} non-selected images to: {}", moved, remove_dir.display());
}

fn divide_by_quality(dir: &Path, images: &[ImageInfo]) {
    use std::collections::BTreeMap;
    let mut bin_counts: BTreeMap<u32, usize> = BTreeMap::new();
    for img in images {
        let bin = (img.quality * 100.0).round() as u32;
        let bin_dir = dir.join(bin.to_string());
        fs::create_dir_all(&bin_dir).expect("Failed to create bin directory");
        let dest = bin_dir.join(&img.file_name);
        if let Err(_) = fs::rename(&img.file_path, &dest) {
            // Cross-device fallback: copy then delete
            fs::copy(&img.file_path, &dest)
                .expect("Failed to copy file to bin directory");
            fs::remove_file(&img.file_path)
                .expect("Failed to remove original after copy");
        }
        *bin_counts.entry(bin).or_insert(0) += 1;
    }
    for (bin, count) in &bin_counts {
        println!("  {}/  <- {} file(s)", bin, count);
    }
    println!("Divided {} files into {} bin(s).", images.len(), bin_counts.len());
}

fn config_fingerprint(config: &AppConfig, crop: Option<f64>, pattern: Option<&[RegisteredStar]>) -> String {
    // Include every input which affects detection/matching.  Pattern JSON is
    // represented by its decoded values, making insignificant JSON formatting
    // changes harmless while changing a star invalidates the cache.
    let pattern_json = pattern.map(|p| p.iter().map(|s| (s.pos.x, s.pos.y, s.magnitude, s.use_in_quality)).collect::<Vec<_>>());
    serde_json::to_string(&(config, crop, pattern_json)).unwrap_or_default()
}

fn image_to_metrics(image: &ImageInfo) -> metrics::MetricValues {
    metrics::MetricValues {
        quality: image.quality, fwhm: image.fwhm, star_count: image.star_count,
        brightest_star_adu: image.brightest_star_adu,
        background_raw_adu: image.background_raw_adu,
        background_corrected_adu: image.background_corrected_adu,
        quality_star_pattern: image.quality_image,
        star_brightness_adu: image.star_brightness_adu,
        snr: image.snr,
        star_pattern_found: image.constellation_found == Some(true),
        matched_star_count: image.matched_star_count,
    }
}

fn record_to_image(record: &metrics::CacheRecord, dir: &Path) -> ImageInfo {
    let m = &record.metrics;
    ImageInfo { file_name: record.file_name.clone(), file_path: dir.join(&record.file_name),
        modified_ns: record.modified_ns, size: record.size, quality: m.quality, fwhm: m.fwhm,
        quality_image: m.quality_star_pattern, star_count: m.star_count,
        constellation_found: if m.star_pattern_found { Some(true) } else { Some(false) },
        matched_star_count: m.matched_star_count, star_brightness_adu: m.star_brightness_adu,
        snr: m.snr, brightest_star_adu: m.brightest_star_adu,
        background_raw_adu: m.background_raw_adu, background_corrected_adu: m.background_corrected_adu,
        brightest_star_photons: 0.0, star5_photons: 0.0 }
}

fn build_rules(args: &Args, has_pattern: bool) -> Result<Vec<metrics::FilterRule>, String> {
    use metrics::Metric::*;
    let mut rules = Vec::new();
    let mut add = |metric: metrics::Metric, relative: Option<f64>, absolute: Option<f64>| -> Result<(), String> {
        if relative.is_some() && absolute.is_some() { return Err(format!("both relative and absolute thresholds supplied for {}", metric.key())); }
        if let Some(v) = relative { rules.push(metrics::FilterRule::relative(metric, v)?); }
        if let Some(v) = absolute { rules.push(metrics::FilterRule::absolute(metric, v)?); }
        Ok(())
    };
    add(Snr, args.snr, args.snr_absolute)?;
    add(QualityStarPattern, args.quality_star_pattern, args.quality_star_pattern_absolute)?;
    add(Background, args.background, args.background_absolute)?;
    add(StarBrightness, args.star_brightness, args.star_brightness_absolute)?;
    add(Quality, args.quality, args.quality_absolute)?;
    add(Fwhm, args.fwhm, args.fwhm_absolute)?;
    if args.filter && rules.is_empty() {
        rules.push(metrics::FilterRule::relative(QualityStarPattern, 0.83)?);
        rules.push(metrics::FilterRule::relative(Snr, 0.707)?);
    }
    if rules.iter().any(|r| r.metric.requires_pattern()) && !has_pattern {
        return Err("filter requires a star pattern, but no pattern is available".into());
    }
    Ok(rules)
}

fn print_statistics(records: &[metrics::CacheRecord]) {
    println!("\n--- METRIC STATISTICS ({} records) ---", records.len());
    for metric in [metrics::Metric::Quality, metrics::Metric::Fwhm, metrics::Metric::QualityStarPattern, metrics::Metric::Background, metrics::Metric::StarBrightness, metrics::Metric::Snr] {
        let mut values: Vec<f64> = records.iter().filter_map(|r| metric.value(&r.metrics)).filter(|v| v.is_finite()).collect();
        if values.is_empty() { continue; }
        let count = values.len(); let mean = values.iter().sum::<f64>() / count as f64;
        values.sort_by(|a,b| a.total_cmp(b));
        println!("{:<24} count={:<4} median={:.4} mean={:.4} min={:.4} max={:.4}", metric.key(), count, metrics::median(&mut values.clone()).unwrap(), mean, values[0], values[count-1]);
    }
    if records.iter().any(|r| !r.metrics.star_pattern_found) { println!("Warning: images without a matched pattern will fail pattern-dependent filters."); }
}

fn draw_metric_chart(dir: &Path, records: &[metrics::CacheRecord], metric: metrics::Metric) -> Result<(), String> {
    use image::{ImageBuffer, Rgb};
    let mut ordered: Vec<&metrics::CacheRecord> = records.iter().collect();
    ordered.sort_by(|a, b| a.modified_ns.cmp(&b.modified_ns).then_with(|| a.file_name.cmp(&b.file_name)));
    let points: Vec<f64> = ordered.iter().filter_map(|r| metric.value(&r.metrics)).collect();
    if points.is_empty() { return Ok(()); }
    let width = 900u32; let height = 420u32;
    let mut out: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::from_pixel(width, height, Rgb([255u8,255u8,255u8]));
    let min = points.iter().copied().fold(f64::INFINITY, f64::min);
    let max = points.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let range = (max - min).max(1e-12);
    for (i, value) in points.iter().enumerate() {
        let x = 20 + (i as u32 * (width - 40) / points.len().max(1) as u32);
        let y = height - 20 - (((*value - min) / range) * (height - 40) as f64) as u32;
        for dx in 0..4 { for dy in 0..4 { if x+dx < width && y+dy < height { out.put_pixel(x+dx,y+dy,Rgb([20,70,180])); } } }
    }
    out.save(dir.join(format!("metrics_{}.png", metric.key()))).map_err(|e| e.to_string())
}

fn filter_files(dir: &Path, records: &[metrics::CacheRecord], rules: &[metrics::FilterRule]) -> Result<(), String> {
    let medians = metrics::medians(records);
    for rule in rules {
        let threshold = rule.threshold(medians.get(&rule.metric).copied())?;
        let median = medians.get(&rule.metric).copied().unwrap_or(0.0);
        println!("{}: median={:.4}, threshold={:.4}, direction {}", rule.metric.key(), median, threshold, if rule.metric.higher_is_better() { ">=" } else { "<=" });
    }
    let rejected: Vec<_> = records.iter().filter(|r| !metrics::passes_all_filters(r, rules, &medians).unwrap_or(false)).collect();
    let destination = metrics::unique_removed_folder(dir, rules, metrics::current_timestamp());
    fs::create_dir_all(&destination).map_err(|e| e.to_string())?;
    for record in &rejected {
        let source = dir.join(&record.file_name); let target = destination.join(&record.file_name);
        fs::rename(&source, &target).or_else(|_| { fs::copy(&source, &target).map(|_| ()).and_then(|_| fs::remove_file(&source)) }).map_err(|e| e.to_string())?;
    }
    println!("Kept {}, moved {} files to {}", records.len() - rejected.len(), rejected.len(), destination.display());
    Ok(())
}

fn process_directory(dir: &Path, args: &Args, registered_stars: Option<&Vec<RegisteredStar>>, config: &AppConfig) {
    let all_entries = collect_fits_files(dir);
    if all_entries.is_empty() { eprintln!("No FITS files found."); return; }
    let fingerprint = config_fingerprint(config, args.crop, registered_stars.map(|s| s.as_slice()));
    let cache_path = dir.join("metrics_cache.json");
    let mut cache = metrics::MetricsCache::load(&cache_path, &fingerprint).unwrap_or_else(|_| metrics::MetricsCache::empty(&fingerprint));
    let entries: Vec<_> = if let Some(n) = args.check_count {
        let mut sorted: Vec<_> = all_entries.iter().map(|e| e.path()).collect(); sorted.sort_by(|a,b| fs::metadata(a).and_then(|m|m.modified()).ok().cmp(&fs::metadata(b).and_then(|m|m.modified()).ok()).then_with(|| a.file_name().cmp(&b.file_name())));
        sorted.into_iter().rev().take(n).filter_map(|p| fs::read_dir(dir).ok()?.find_map(|e| e.ok().filter(|e| e.path() == p))).collect()
    } else { all_entries.into_iter().collect() };
    let stale_paths: Vec<PathBuf> = entries.iter().filter(|e| !cache.record_is_current(&e.path(), &fingerprint)).map(|e| e.path()).collect();
    let stale: Vec<_> = stale_paths.iter().filter_map(|path| fs::read_dir(dir).ok()?.find_map(|e| e.ok().filter(|e| e.path() == *path))).collect();
    let analyzed = load_images(&stale, args.crop, config, registered_stars);
    for image in &analyzed {
        if let Ok(record) = metrics::record_for(&image.file_path, &fingerprint, image_to_metrics(image)) { cache.upsert(record); }
    }
    if let Err(e) = cache.save(&cache_path) { eprintln!("Warning: cannot save cache: {}", e); }
    let records: Vec<_> = if args.check_count.is_some() { entries.iter().filter_map(|e| cache.records.iter().find(|r| r.file_name == e.file_name().to_string_lossy())).cloned().collect() } else { cache.records.iter().filter(|r| dir.join(&r.file_name).is_file()).cloned().collect() };
    let images: Vec<_> = records.iter().map(|r| record_to_image(r, dir)).collect();
    write_quality_map(dir, &images, None, config);
    print_statistics(&records);
    if args.check_count.is_none() { println!("Suggested filters: quality_star_pattern >= 0.83 × median; snr >= 0.707 × median"); }
    if args.check_count.is_none() && args.filter == false && build_rules(args, registered_stars.is_some()).map(|r| !r.is_empty()).unwrap_or(false) { eprintln!("Warning: thresholds supplied without --filter; no files moved."); }
    if args.check_count.is_none() && args.filter {
        match build_rules(args, registered_stars.is_some()).and_then(|rules| filter_files(dir, &records, &rules)) { Ok(()) => {}, Err(e) => eprintln!("Error: {}", e) }
    }
    if args.check_count.is_none() { for metric in [metrics::Metric::Quality, metrics::Metric::Fwhm, metrics::Metric::QualityStarPattern, metrics::Metric::Background, metrics::Metric::StarBrightness, metrics::Metric::Snr] { let _ = draw_metric_chart(dir, &records, metric); } }
}

fn check_seeing(args: &Args, registered_stars: &Option<Vec<RegisteredStar>>, config: &AppConfig, path: &Path) {
    let entries = collect_fits_files(path);
    if entries.is_empty() {
        eprintln!("No FITS files found in directory for seeing check.");
        std::process::exit(1);
    }

    let sample_size = args.check_seeing.expect("seeing sample size is required");
    if sample_size == 0 {
        eprintln!("Seeing sample size must be greater than zero.");
        std::process::exit(1);
    }

    // Use creation time so repeated calls such as `-c 3` analyze the latest
    // batch of frames rather than an arbitrary subset of the directory.
    let mut entries = entries;
    entries.sort_by(|a, b| {
        let a_created = a.metadata().and_then(|metadata| metadata.modified()).ok();
        let b_created = b.metadata().and_then(|metadata| metadata.modified()).ok();
        a_created
            .cmp(&b_created)
            .then_with(|| a.file_name().cmp(&b.file_name()))
    });
    let start = entries.len().saturating_sub(sample_size);
    let sampled_entries = &entries[start..];

    let images = load_images(sampled_entries, args.crop, &config, registered_stars.as_ref());
    if images.is_empty() {
        eprintln!("No images loaded for seeing check.");
        std::process::exit(1);
    }
    let mut qualities: Vec<f64> = images.iter().map(|i| i.quality).collect();
    qualities.sort_by(|a, b| a.total_cmp(b));
    let mean = qualities.iter().sum::<f64>() / qualities.len() as f64;
    let stddev = if qualities.len() > 1 {
        let mean = mean;
        (qualities.iter().map(|q| (q - mean).powi(2)).sum::<f64>() / (qualities.len() as f64 - 1.0)).sqrt()
    } else { 0.0 };
    let median = percentile(&qualities, 50);
    let pct25 = percentile(&qualities, 25);
    let pct75 = percentile(&qualities, 75);
    println!("\n--- SEEING METRICS ({} images) ---", qualities.len());
    println!("Median FWHM (SEEING): {:.4}", fwhm_from_quality(median));
    println!("FWHM VARY FROM: {:.4} - {:.4}", fwhm_from_quality(mean + stddev), fwhm_from_quality(mean - stddev));
    println!("Median QUALITY: {:.4}%", median * 100.0);
    println!("Mean QUALITY: {:.4}%", mean * 100.0);
    println!("Stddev QUALITY: {:.4}%", stddev * 100.0);
    println!("25th percentile: {:.4}%", pct25 * 100.0);
    println!("75th percentile: {:.4}%", pct75 * 100.0);
    println!("-----------------------------------\n");
}

fn make_star_pattern(folder: &Path, config: &AppConfig) -> Result<(), Box<dyn std::error::Error>> {
    let entries = collect_fits_files(folder);
    if entries.is_empty() { return Err("folder contains no FITS files".into()); }
    let sample_paths = star_pattern::select_sample_paths(entries.iter().map(|e| e.path()).collect(), 20, 0);
    let mut samples = Vec::new();
    let mut first: Option<AstroImage> = None;
    for path in sample_paths {
        let img = AstroImage::load(&path, None, config)?;
        if first.is_none() { first = Some(img); } else {
            for s in img.stars() { samples.push(star_pattern::PatternStarSample { x:s.pos.x, y:s.pos.y, magnitude:s.magnitude, brightest_pixel_part:s.brightest_pixel_part }); }
        }
    }
    let image = first.ok_or("unable to load FITS samples")?;
    let rec = star_pattern::recommend_stars(&samples, image.width(), image.height());
    if rec.len() != 3 { return Err("input data does not contain a sensible three-star pattern".into()); }
    let mut pattern = Vec::new();
    for &idx in &rec {
        let s = &samples[idx];
        pattern.push(star_pattern::StarPatternEntry { x:s.x, y:s.y, magnitude:s.magnitude, use_in_quality:true, median_brightness:Some(s.magnitude), median_brightest_pixel_part:Some(s.brightest_pixel_part) });
    }
    image.save_stars_jpg(folder.join("star_pattern_candidates.jpg"))?;
    println!("Recommended stars: {:?}", rec.iter().map(|i| i + 1).collect::<Vec<_>>());
    println!("Accept recommendation? [Y/n]");
    let mut answer = String::new(); std::io::stdin().read_line(&mut answer)?;
    if answer.trim().eq_ignore_ascii_case("n") { return Err("manual selection is not available without --star-pattern-numbers".into()); }
    let output = folder.join("stars_pattern.json");
    std::fs::write(&output, serde_json::to_string_pretty(&pattern)?)?;
    println!("Pattern written to {}", output.display());
    Ok(())
}

fn main() {
    let args = Args::parse();

    let registered_stars: Option<Vec<RegisteredStar>> = match &args.star_pattern {
        Some(p) => Some(load_stars_from_json(p).unwrap_or_else(|e| {
            eprintln!("Error loading star pattern '{}': {}", p, e);
            std::process::exit(1);
        })),
        None => {
            let default_stars = star_pattern::default_pattern_path(Path::new(&args.path));
            if let Some(default_stars) = default_stars {
                match load_stars_from_json(default_stars) {
                    Ok(stars) => {
                        println!("Using stars.json from current directory.");
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
    if let Err(e) = config.validate() { eprintln!("Error: {}", e); std::process::exit(2); }
    if let Some(folder) = &args.make_star_pattern {
        if let Err(e) = make_star_pattern(Path::new(folder), &config) { eprintln!("Error creating star pattern: {}", e); std::process::exit(1); }
        return;
    }

    let path = Path::new(&args.path);
    if path.is_file() {
        process_single_file(path, args.crop, args.save_stars, &config, registered_stars.as_ref());
    } else if path.is_dir() {
        if args.check_seeing.is_some() {
            check_seeing(&args, &registered_stars, &config, path);
        } else if args.save_stars {
            let entries = collect_fits_files(path);
            for entry in entries {
                process_single_file(&entry.path(), args.crop, true, &config, registered_stars.as_ref());
            }
        } else {
            process_directory(path, &args, registered_stars.as_ref(), &config);
        }
    } else {
        eprintln!("Error: path does not exist: {}", args.path);
        std::process::exit(1);
    }
}



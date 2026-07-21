use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use serde::Deserialize;
use crate::astro_image::{fwhm_from_quality, AstroImage};
use crate::constellation::{Constellation, RegisteredStar, load_stars_from_json};
use crate::helpers::median;

pub mod astro_image;
pub mod star;
pub mod helpers;
pub mod constellation;


pub struct AppConfig {
    gain_to_adu: HashMap<u32, f64>,
    min_photons_to_detect_star: i32,
    min_central_photons_to_detect_star: i32,
    psf_size: usize,
    min_photons_quality: f64,
    rolling_avg_window: usize,
    log_quality_window_t: f64,
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
        })
    }
}

impl AppConfig {
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

    /// Copy the top N best images to a 'selected' subfolder (folder mode only)
    #[arg(long, short, value_name = "FRACTION")]
    take: Option<f64>,

    /// Copy the top N best images based on quality_image score to a named subfolder (folder mode only)
    #[arg(long, value_name = "FRACTION")]
    take_quality: Option<f64>,

    /// Move non-selected images to a 'remove' subfolder instead of deleting them (folder mode only)
    #[arg(long, short)]
    remove: bool,

    /// Only search for stars in the central fraction of the image (e.g. 0.3 = central 30% width and height)
    #[arg(long, value_name = "FRACTION")]
    crop: Option<f64>,

    /// Save annotated star image to a JPG file
    #[arg(long, short)]
    save_stars: bool,

    /// Path to a JSON file with reference stars for constellation-based quality filtering
    #[arg(long, value_name = "FILE")]
    star_pattern: Option<String>,

    /// Keep images with low stars
    #[arg(long, short)]
    keep_low: bool,

    /// Check seeing using the N most recently created FITS images
    #[arg(long, short = 'c', value_name = "N")]
    check_seeing: Option<usize>,

    /// Move all files into subdirectories named by their rounded quality percentage (folder mode only)
    #[arg(long)]
    divide: bool,
}

struct ImageInfo {
    file_name: String,
    file_path: PathBuf,
    quality: f64,
    fwhm: f64,
    quality_image: Option<f64>,
    star_count: usize,
    constellation_found: Option<bool>,
    brightest_star_photons: f64,
    star5_photons: f64,
}



fn apply_constellation_quality(img: &mut AstroImage, registered_stars: &[RegisteredStar], config: &AppConfig, label: &str) -> bool {
    let constellation = Constellation::find_in_image(registered_stars.to_vec(), img);
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
                if !was_images_in_window_set {
                    images_in_window = (config.log_quality_window_t / img.exp_t()).ceil() as i32;
                    was_images_in_window_set = true
                }
                let constellation_found = registered_stars.map(|stars| {
                    apply_constellation_quality(&mut img, stars, config, &file_name)
                });
                if img.quality().is_finite(){
                    quality_sum += img.quality();
                    fwhm_sum += img.fwhm();
                    quality_count += 1;
                    last_images_window.push(img.quality());
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
                    fwhm: img.fwhm(),
                    quality: img.quality(),
                    quality_image: img.quality_image(),
                    star_count: img.star_count(),
                    file_name,
                    file_path,
                    constellation_found,
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

fn process_directory(dir: &Path, args: &Args, registered_stars: Option<&Vec<RegisteredStar>>, config: &AppConfig) {
    if args.divide && (args.take.is_some() || args.take_quality.is_some()) {
        eprintln!("Error: --divide cannot be combined with --take or --take-quality.");
        return;
    }

    let entries = collect_fits_files(dir);
    let images = load_images(&entries, args.crop, config, registered_stars);

    let threshold_info = if !images.is_empty() { Some(compute_star_threshold(&images)) } else { None };
    write_quality_map(dir, &images, threshold_info.map(|t| t.1), config);

    let use_constellation = registered_stars.is_some();
    let mut all_selected: HashSet<&str> = HashSet::new();

    if let Some(take_pct) = args.take {
        if images.is_empty() {
            eprintln!("No images loaded.");
            return;
        }
        let (median_stars, low_star_threshold) = threshold_info.unwrap();
        let pct_int = (take_pct * 100.0).round() as u32;
        let folder_name = format!("selected_percent_{}", pct_int);
        let selected = select_best_by_percent(&images, take_pct, median_stars, low_star_threshold, use_constellation, args.keep_low);
        copy_to_named_folder(dir, &images, &selected, &folder_name);
        all_selected.extend(selected);
    }

    if let Some(take_quality) = args.take_quality {
        if images.is_empty() {
            eprintln!("No images loaded.");
            return;
        }
        let pct_int = (take_quality * 100.0).round() as u32;
        let folder_name = format!("selected_quality_{}", pct_int);
        let selected = select_best_by_quality(&images, take_quality, use_constellation);
        copy_to_named_folder(dir, &images, &selected, &folder_name);
        all_selected.extend(selected);
    }

    if args.remove {
        if args.take.is_none() && args.take_quality.is_none() {
            eprintln!("Warning: --remove has no effect without --take or --take-quality.");
        } else {
            remove_original_images(dir, &images, &all_selected);
        }
    }

    if args.divide {
        divide_by_quality(dir, &images);
    }
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

fn main() {
    let args = Args::parse();

    let registered_stars: Option<Vec<RegisteredStar>> = match &args.star_pattern {
        Some(p) => Some(load_stars_from_json(p).unwrap_or_else(|e| {
            eprintln!("Error loading star pattern '{}': {}", p, e);
            std::process::exit(1);
        })),
        None => {
            let default_stars = Path::new("stars.json");
            if default_stars.exists() {
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



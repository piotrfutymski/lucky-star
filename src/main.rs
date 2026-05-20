use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use serde::Deserialize;
use crate::astro_image::AstroImage;
use crate::constellation::{Constellation, RegisteredStar, load_stars_from_json};
use crate::helpers::median;

pub mod astro_image;
pub mod star;
pub mod helpers;
pub mod constellation;


#[derive(Deserialize)]
pub struct AppConfig {
    gain_to_adu: HashMap<u32, f64>,
    min_photons_to_detect_star: i32,
    min_central_photons_to_detect_star: i32,
    psf_size: usize,
    min_photons_quality: f64,
    rolling_avg_window: usize,
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
            min_photons_to_detect_star: 150,
            min_central_photons_to_detect_star: 12,
            psf_size: 13,
            min_photons_quality: 200.0,
            rolling_avg_window: 100,
        }
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
    /// Path to a FITS file or a directory containing FITS files
    #[arg(default_value = ".")]
    path: String,

    /// Copy the top N best images to a 'selected' subfolder (folder mode only)
    #[arg(long, value_name = "FRACTION")]
    take: Option<f64>,

    /// Copy the top N best images based on quality_image score to a named subfolder (folder mode only)
    #[arg(long, value_name = "FRACTION")]
    take_quality: Option<f64>,

    /// Move non-selected images to a 'remove' subfolder instead of deleting them (folder mode only)
    #[arg(long)]
    remove: bool,

    /// Only search for stars in the central fraction of the image (e.g. 0.3 = central 30% width and height)
    #[arg(long, value_name = "FRACTION")]
    crop: Option<f64>,

    /// Save annotated star image to a JPG file (single file mode only)
    #[arg(long, short)]
    save_stars: bool,

    /// Path to a JSON file with reference stars for constellation-based quality filtering
    #[arg(long, value_name = "FILE")]
    star_pattern: Option<String>
}

struct ImageInfo {
    file_name: String,
    file_path: PathBuf,
    quality: f64,
    quality_image: Option<f64>,
    star_count: usize,
    constellation_found: Option<bool>,
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

fn load_images(entries: Vec<fs::DirEntry>, crop: Option<f64>, config: &AppConfig, registered_stars: Option<&Vec<RegisteredStar>>) -> Vec<ImageInfo> {
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
    let mut quality_image_sum = 0.0f64;
    let mut quality_image_count = 0usize;

    for entry in entries {
        let file_path = entry.path();
        let file_name = entry.file_name().to_string_lossy().to_string();
        match AstroImage::load(&file_path, crop, config) {
            Ok(mut img) => {
                let constellation_found = registered_stars.map(|stars| {
                    apply_constellation_quality(&mut img, stars, config, &file_name)
                });
                quality_sum += img.quality();
                if let Some(qi) = img.quality_image() {
                    quality_image_sum += qi;
                    quality_image_count += 1;
                }
                let n = images.len() + 1;
                let avg_q = quality_sum / n as f64;
                let msg = if quality_image_count > 0 {
                    let avg_qi = quality_image_sum / quality_image_count as f64;
                    format!(
                        "avg quality: {:.2}%  avg quality_image: {:.2}%",
                        avg_q * 100.0,
                        avg_qi * 100.0
                    )
                } else {
                    format!("avg quality: {:.2}%", avg_q * 100.0)
                };
                pb.set_message(msg);
                pb.inc(1);
                images.push(ImageInfo {
                    quality: img.quality(),
                    quality_image: img.quality_image(),
                    star_count: img.star_count(),
                    file_name,
                    file_path,
                    constellation_found,
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

    writeln!(map_file, "# Percentiles").expect("Failed to write");
    writeln!(map_file, "# pct\tquality\tquality_image").expect("Failed to write");
    for p in (0..=100).step_by(10) {
        let q = percentile(&quality_vals, p);
        let qi = if quality_image_vals.is_empty() {
            String::new()
        } else {
            format!("{:.6}", percentile(&quality_image_vals, p))
        };
        if p == 50 {
            println!("MEDIAN QUALITY FOR SEQUENCE: {:.4} %", q * 100.0)
        }
        writeln!(map_file, "# {:>3}%\t{:.6}\t{}", p, q, qi).expect("Failed to write percentile");
    }
    writeln!(map_file, "#").expect("Failed to write");
    writeln!(map_file, "filename\tquality\tquality_image\tstars\tnote").expect("Failed to write header");

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
            .map(|q| format!("{:.6}", q))
            .unwrap_or_default();
        writeln!(
            map_file,
            "{}\t{:.6}\t{}\t{}\t{}",
            img.file_name, img.quality, quality_image_str, img.star_count, note
        )
        .expect("Failed to write quality map");
    }

    // Trend section: rolling averages sorted by filename (chronological proxy)
    let window = config.rolling_avg_window;
    if window > 0 && images.len() > window {
        let mut by_filename: Vec<&ImageInfo> = images.iter().collect();
        by_filename.sort_by(|a, b| a.file_name.cmp(&b.file_name));

        writeln!(map_file, "\n# Trend (rolling average, window = {})", window)
            .expect("Failed to write trend header");
        writeln!(map_file, "# images\tquality\tquality_image")
            .expect("Failed to write trend header");

        let mut start = 0;
        while start < by_filename.len() {
            let end = (start + window).min(by_filename.len());
            let slice = &by_filename[start..end];
            let avg_q = slice.iter().map(|i| i.quality).sum::<f64>() / slice.len() as f64;
            let qi_vals: Vec<f64> = slice.iter().filter_map(|i| i.quality_image).collect();
            let avg_qi_str = if qi_vals.is_empty() {
                String::new()
            } else {
                format!("{:.6}", qi_vals.iter().sum::<f64>() / qi_vals.len() as f64)
            };
            writeln!(map_file, "# {}-{}\t{:.6}\t{}", start + 1, end, avg_q, avg_qi_str)
                .expect("Failed to write trend row");
            start += window;
        }
    }

    println!("\nQuality map written to: {}", map_path.display());
}

fn compute_star_threshold(images: &[ImageInfo]) -> (usize, usize) {
    let mut sorted_stars: Vec<usize> = images.iter().map(|i| i.star_count).collect();
    sorted_stars.sort_unstable();
    let median = median(&sorted_stars).unwrap_or_default();
    let threshold = (median as f64 * 0.7) as usize;
    (median, threshold)
}
fn select_best_by_quality(images: &[ImageInfo], take_pct: f64, median_stars: usize, low_star_threshold: usize, use_constellation: bool) -> HashSet<&str> {
    let take_pct = take_pct.clamp(0.0, 1.0);
    let total = images.len();
    let count_to_take = ((total as f64 * take_pct).ceil() as usize)
        .max(1)
        .min(total);

    let mut eligible: Vec<&ImageInfo> = if use_constellation {
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

fn select_best_by_quality_image(images: &[ImageInfo], take_pct: f64, use_constellation: bool) -> HashSet<&str> {
    let take_pct = take_pct.clamp(0.0, 1.0);
    let total = images.len();
    let count_to_take = ((total as f64 * take_pct).ceil() as usize)
        .max(1)
        .min(total);

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
        .take(count_to_take)
        .map(|i| i.file_name.as_str())
        .collect();

    println!(
        "Quality-image selection: {}/{} eligible, copying top {} by quality_image.",
        eligible.len(), total, selected.len()
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

fn process_directory(dir: &Path, args: &Args, registered_stars: Option<&Vec<RegisteredStar>>, config: &AppConfig) {
    let entries = collect_fits_files(dir);
    let images = load_images(entries, args.crop, config, registered_stars);

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
        let selected = select_best_by_quality(&images, take_pct, median_stars, low_star_threshold, use_constellation);
        copy_to_named_folder(dir, &images, &selected, &folder_name);
        all_selected.extend(selected);
    }

    if let Some(take_quality_pct) = args.take_quality {
        if images.is_empty() {
            eprintln!("No images loaded.");
            return;
        }
        let pct_int = (take_quality_pct * 100.0).round() as u32;
        let folder_name = format!("selected_quality_{}", pct_int);
        let selected = select_best_by_quality_image(&images, take_quality_pct, use_constellation);
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
        process_directory(path, &args, registered_stars.as_ref(), &config);
    } else {
        eprintln!("Error: path does not exist: {}", args.path);
        std::process::exit(1);
    }
}

use std::collections::HashSet;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use clap::Parser;
use crate::astro_image::AstroImage;
use crate::helpers::median;

pub mod astro_image;
pub mod star;
pub mod constants;
pub mod helpers;

#[derive(Parser)]
#[command(name = "lucky-star", about = "Astronomical image quality analyzer")]
struct Args {
    /// Path to a FITS file or a directory containing FITS files
    path: String,

    /// Copy the top N best images to a 'selected' subfolder (folder mode only)
    #[arg(long, value_name = "FRACTION")]
    take: Option<f64>,

    /// Remove images not selected by --take (folder mode only)
    #[arg(long)]
    remove: bool,

    /// Only search for stars in the central fraction of the image (e.g. 0.3 = central 30% width and height)
    #[arg(long, value_name = "FRACTION")]
    crop: Option<f64>,

    /// Save annotated star image to a JPG file (single file mode only)
    #[arg(long)]
    save_stars: bool,

    /// Min flux in photons to account for quality
    #[arg(long, default_value_t = 200.0)]
    min_photons_quality: f64,

    /// PSF window size in pixels used for star detection and flux measurement
    #[arg(long, default_value_t = 13)]
    psf_size: usize
}

struct ImageInfo {
    file_name: String,
    file_path: PathBuf,
    quality: f64,
    star_count: usize,
}



fn process_single_file(path: &Path, crop: Option<f64>, save_stars: bool, min_photons_quality: f64, psf_size: usize) {
    let img = AstroImage::load(path, crop, min_photons_quality, psf_size).unwrap();
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

fn load_images(entries: Vec<fs::DirEntry>, crop: Option<f64>, min_photons_quality: f64, psf_size: usize) -> Vec<ImageInfo> {
    let mut images = Vec::new();
    for entry in entries {
        let file_path = entry.path();
        let file_name = entry.file_name().to_string_lossy().to_string();
        match AstroImage::load(&file_path, crop, min_photons_quality, psf_size) {
            Ok(img) => {
                print!("{}", img.brief_summary());
                images.push(ImageInfo {
                    quality: img.quality(),
                    star_count: img.star_count(),
                    file_name,
                    file_path,
                });
            }
            Err(e) => eprintln!("Error loading {}: {}", file_name, e),
        }
    }
    images
}

fn write_quality_map(dir: &Path, images: &[ImageInfo], low_star_threshold: Option<usize>) {
    let map_path = dir.join("quality_map.txt");
    let mut map_file = fs::File::create(&map_path).expect("Failed to create quality_map.txt");
    writeln!(map_file, "filename\tquality\tstars\tnote").expect("Failed to write header");

    let mut sorted: Vec<&ImageInfo> = images.iter().collect();
    sorted.sort_by(|a, b| b.quality.total_cmp(&a.quality));
    for img in &sorted {
        let note = match low_star_threshold {
            Some(t) if img.star_count < t => "low_stars",
            _ => "",
        };
        writeln!(map_file, "{}\t{:.6}\t{}\t{}", img.file_name, img.quality, img.star_count, note)
            .expect("Failed to write quality map");
    }
    println!("\nQuality map written to: {}", map_path.display());
}

fn compute_star_threshold(images: &[ImageInfo]) -> (usize, usize) {
    let mut sorted_stars: Vec<usize> = images.iter().map(|i| i.star_count).collect();
    sorted_stars.sort_unstable();
    let median = median(&sorted_stars).unwrap_or_default();
    let threshold = (median as f64 * 0.7) as usize; //TODO CHANGE LOGIC IN FUTURE
    (median, threshold)
}

fn select_best_images(images: &[ImageInfo], take_pct: f64, median_stars: usize, low_star_threshold: usize) -> HashSet<&str> {
    let take_pct = take_pct.clamp(0.0, 100.0);
    let total = images.len();
    let count_to_take = ((total as f64 * take_pct / 100.0).ceil() as usize)
        .max(1)
        .min(total);

    let mut eligible: Vec<&ImageInfo> = images.iter()
        .filter(|i| i.star_count >= low_star_threshold)
        .collect();
    eligible.sort_by(|a, b| b.quality.total_cmp(&a.quality));

    let selected: HashSet<&str> = eligible.iter()
        .take(count_to_take)
        .map(|i| i.file_name.as_str())
        .collect();

    println!(
        "\nSelection: {}/{} eligible (median {} stars, min threshold {}), copying top {} by quality.",
        eligible.len(), total, median_stars, low_star_threshold, selected.len()
    );

    selected
}

fn copy_selected_images(dir: &Path, images: &[ImageInfo], selected: &HashSet<&str>) {
    let selected_dir = dir.join("selected");
    fs::create_dir_all(&selected_dir).expect("Failed to create 'selected' directory");
    for img in images {
        if selected.contains(img.file_name.as_str()) {
            let dest = selected_dir.join(&img.file_name);
            fs::copy(&img.file_path, &dest).expect("Failed to copy file");
        }
    }
    println!("Copied {} images to: {}", selected.len(), selected_dir.display());
}

fn remove_original_images(images: &[ImageInfo]) {
    let mut removed = 0usize;
    for img in images {
        fs::remove_file(&img.file_path).expect("Failed to remove file");
        removed += 1;
    }
    println!("Removed {} non-selected images.", removed);
}

fn process_directory(dir: &Path, args: &Args) {
    let entries = collect_fits_files(dir);
    let images = load_images(entries, args.crop, args.min_photons_quality, args.psf_size);

    let threshold_info = if !images.is_empty() { Some(compute_star_threshold(&images)) } else { None };
    write_quality_map(dir, &images, threshold_info.map(|info| info.1));

    if let Some(take_pct) = args.take {
        if images.is_empty() {
            eprintln!("No images loaded.");
            return;
        }
        let (median_stars, low_star_threshold) = threshold_info.unwrap();
        let selected = select_best_images(&images, take_pct, median_stars, low_star_threshold);
        copy_selected_images(dir, &images, &selected);
        if args.remove {
            remove_original_images(&images);
        }
    } else if args.remove {
        eprintln!("Warning: --remove has no effect without --take.");
    }
}

fn main() {
    let args = Args::parse();
    let path = Path::new(&args.path);

    if path.is_file() {
        process_single_file(path, args.crop, args.save_stars, args.min_photons_quality, args.psf_size);
    } else if path.is_dir() {
        process_directory(path, &args);
    } else {
        eprintln!("Error: path does not exist: {}", args.path);
        std::process::exit(1);
    }
}

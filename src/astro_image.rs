use crate::AppConfig;
use crate::star::Star;
use ndarray::{Array, Array1, IxDyn};
use rand::RngExt;
use rustronomy_fits::{Extension, Fits, Header};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::fs;
use std::io::Write;
use std::path::Path;
use vector2d::Vector2D;

pub struct AstroImage {
    width: u32,
    height: u32,
    gain: u32,
    adu_e: f64,
    exp_t: f64,

    background_level_adu: u16,
    sigma_adu: u16,

    psf_size: usize,
    detected_stars: Vec<Star>,
    data: Array<i16, IxDyn>,
    fwhm: f64,
    quality: f64,
    quality_image: Option<f64>,
    quality_star_indices: Option<Vec<usize>>,
}

impl AstroImage {
    pub fn load(
        file: impl AsRef<Path>,
        crop: Option<f64>,
        config: &AppConfig,
    ) -> Result<AstroImage, Box<dyn Error>> {
        let fits = Fits::open(file.as_ref())?;

        let error_msg = format!("Can not get data from fit file: {:?}", file.as_ref());
        let header = fits.get_hdu(0).expect(error_msg.as_str()).get_header();
        let mut res = Self::get_data_from_header(header, config)?;
        res.psf_size = config.psf_size;

        let data_array = match fits.get_hdu(0).expect(error_msg.as_str()).get_data() {
            Some(Extension::Image(img)) => img.as_i16_array()?,
            _ => Err(error_msg.clone())?,
        };

        res.detect_stars(data_array, crop, config);
        res.data = data_array.clone();
        Ok(res)
    }

    fn get_data_from_header(
        header: &Header,
        config: &AppConfig,
    ) -> Result<AstroImage, Box<dyn Error>> {
        let width: u32 = header.get_value_as("NAXIS1")?;
        let height: u32 = header.get_value_as("NAXIS2")?;
        let exp_t: f64 = header.get_value_as("EXPOSURE").unwrap_or(0.0);
        let gain: u32 = header.get_value_as("GAIN").unwrap_or(0);
        let adu_e: f64 = *config
            .gain_to_adu
            .get(&gain)
            .unwrap_or_else(|| panic!("No adu_e defined for gain {}", gain));
        let res = Self {
            width,
            height,
            gain,
            adu_e,
            exp_t,
            background_level_adu: 0,
            sigma_adu: 0,
            detected_stars: vec![],
            data: Array::zeros(IxDyn(&[0, 0])),
            quality: 0.0,
            quality_image: None,
            psf_size: 0,
            quality_star_indices: None,
            fwhm: 0.0,
        };
        Ok(res)
    }

    fn detect_stars(&mut self, data: &Array<i16, IxDyn>, crop: Option<f64>, config: &AppConfig) {
        self.extract_global_image_metadata(data);
        self.find_stars(data, crop, config);
        let mut quality_sum = 0.0;
        let mut mag_sum = 0.0;
        self.detected_stars
            .iter()
            .filter(|star| star.magnitude > config.min_photons_quality)
            .for_each(|s| {
                quality_sum += s.magnitude * s.top_4_pixels_part;
                mag_sum += s.magnitude;
            });
        if mag_sum > 0.0 {
            self.quality = quality_sum / mag_sum;
            self.recalculate_fwhm();
        } else {
            self.quality = 0.0;
        }
    }

    fn extract_global_image_metadata(&mut self, data: &Array<i16, IxDyn>) {
        let sample_count: usize = 50000;
        let mut rng = rand::rng();
        let samples = (0..sample_count)
            .map(|_| {
                (
                    (rng.random::<u32>() % self.width) as usize,
                    (rng.random::<u32>() % self.height) as usize,
                )
            })
            .map(|idx| (data[[idx.0, idx.1]].wrapping_sub(i16::MIN) as u16) as f32)
            .collect::<Array1<f32>>();

        let mean = samples.mean().unwrap_or(0.);
        let sigma = samples.std(0.);
        let samples_no_stars = samples
            .iter()
            .filter(|s| **s < mean + 3.0 * sigma)
            .copied()
            .collect::<Array1<f32>>();
        let background_value = samples_no_stars.mean().unwrap_or(0.0);
        let sigma = samples_no_stars.std(0.);

        self.background_level_adu = background_value as u16;
        self.sigma_adu = sigma as u16;
    }

    fn find_stars(&mut self, data: &Array<i16, IxDyn>, crop: Option<f64>, config: &AppConfig) {
        let psf_size = self.psf_size;
        let min_v = self.background_level_adu
            + (self.adu_e * config.min_central_photons_to_detect_star as f64) as u16;
        let max_v = 0.7 * u16::MAX as f64;
        let mut potential_stars = HashMap::new();
        let (i_start, i_end, j_start, j_end) = if let Some(c) = crop {
            let c = c.clamp(0.0, 1.0);
            let to_include = (self.width.max(self.height) as f64 * c) as i32;
            let margin_w = ((self.width as i32 - to_include) / 2).max(psf_size as i32) as usize;
            let margin_h = ((self.height as i32 - to_include) / 2).max(psf_size as i32) as usize;
            let i_end = (self.width as usize - margin_w).min(self.width as usize - psf_size);
            let j_end = (self.height as usize - margin_h).min(self.height as usize - psf_size);
            (margin_w, i_end, margin_h, j_end)
        } else {
            (
                psf_size,
                self.width as usize - psf_size,
                psf_size,
                self.height as usize - psf_size,
            )
        };
        for i in i_start..i_end {
            for j in j_start..j_end {
                let v = data[[i, j]].wrapping_sub(i16::MIN) as u16;
                if v > min_v {
                    potential_stars.insert((i, j), (v, true));
                }
            }
        }
        for key in potential_stars.keys().cloned().collect::<Vec<_>>() {
            let potential = potential_stars.get(&key).cloned().unwrap();
            if potential.1 {
                let mut is_star = true;
                for i in (key.0 - 2)..=(key.0 + 2) {
                    for j in (key.1 - 2)..=(key.1 + 2) {
                        if i == key.0 && j == key.1 {
                            continue;
                        }
                        if let Some(potential_next_to) = potential_stars.get_mut(&(i, j)) {
                            if potential_next_to.0 >= potential.0 {
                                is_star = false;
                            } else {
                                potential_next_to.1 = false;
                            }
                        }
                    }
                }
                if !is_star {
                    potential_stars.get_mut(&key).unwrap().1 = false;
                }
            }
        }
        let mut star_pos = potential_stars
            .iter()
            .filter(|(_k, v)| v.1)
            .map(|(k, _v)| *k)
            .collect::<Vec<_>>();
        let mut to_delete = HashSet::new();
        for (i, pos) in star_pos.iter().enumerate() {
            for (j, compare_pos) in star_pos.iter().enumerate() {
                if j <= i {
                    continue;
                }
                if (compare_pos.0 as i32 - pos.0 as i32).pow(2)
                    + (compare_pos.1 as i32 - pos.1 as i32).pow(2)
                    < (psf_size * psf_size) as i32
                {
                    to_delete.insert(i);
                    to_delete.insert(j);
                    break;
                }
            }
        }
        let mut to_delete = to_delete.into_iter().collect::<Vec<_>>();
        to_delete.sort_by(|a, b| b.cmp(a));
        for index in to_delete {
            star_pos.remove(index);
        }
        let mut stars = star_pos
            .iter()
            .map(|&i| {
                Star::new(
                    Vector2D::new(i.0, i.1),
                    data,
                    self.adu_e,
                    self.background_level_adu,
                    psf_size,
                )
            })
            .filter(|s| {
                !s.ill_defined
                    && s.magnitude > config.min_photons_to_detect_star as f64
                    && s.brightest_pixel_adu < max_v
            })
            .collect::<Vec<_>>();
        if !stars.is_empty() {
            let mut top4_vals: Vec<f64> = stars.iter().map(|s| s.top_4_pixels_part).collect();
            top4_vals.sort_by(|a, b| a.total_cmp(b));
            let n = top4_vals.len();
            let median = if n.is_multiple_of(2) {
                (top4_vals[n / 2 - 1] + top4_vals[n / 2]) / 2.0
            } else {
                top4_vals[n / 2]
            };
            stars.retain(|s| s.top_4_pixels_part >= median / 2.0);
        }
        stars.sort_by(|a, b| {
            b.magnitude
                .partial_cmp(&a.magnitude)
                .unwrap_or(Ordering::Equal)
        });
        self.detected_stars = stars;
    }

    fn recalculate_fwhm(&mut self) {
        self.fwhm = fwhm_from_quality(self.quality)
    }
}

// PUBLIC

pub fn fwhm_from_quality(quality: f64) -> f64 {
    -11.06713 + 26.141212 * (quality * 100.0).powf(-0.186645)
}

impl AstroImage {
    pub fn save_stars_jpg(&self, output_path: impl AsRef<Path>) -> Result<(), Box<dyn Error>> {
        let black_point = self.background_level_adu as f32;
        let white_point = black_point + 20.0 * self.sigma_adu as f32;
        let range = (white_point - black_point).max(1.0);

        let mut img = image::RgbImage::new(self.width, self.height);

        for x in 0..self.width as usize {
            for y in 0..self.height as usize {
                let raw = self.data[[x, y]].wrapping_sub(i16::MIN) as u16 as f32;
                let normalized = ((raw - black_point) / range).clamp(0.0, 1.0);
                let val = (normalized * 255.0) as u8;
                img.put_pixel(x as u32, y as u32, image::Rgb([val, val, val]));
            }
        }

        for (i, star) in self.detected_stars.iter().enumerate() {
            let cx = star.pos.x as i32;
            let cy = star.pos.y as i32;
            imageproc::drawing::draw_hollow_circle_mut(
                &mut img,
                (cx, cy),
                5,
                image::Rgb([255u8, 0u8, 0u8]),
            );
            draw_number_label(&mut img, cx, cy, i + 1);
        }

        img.save(output_path.as_ref())?;
        Ok(())
    }

    pub fn save_stars_md(&self, output_path: impl AsRef<Path>) -> Result<(), Box<dyn Error>> {
        let mut file = fs::File::create(output_path.as_ref())?;
        writeln!(file, "# Detected Stars (by mag)\n")?;
        writeln!(
            file,
            "| # | X | Y | Flux (e\u{207b}) | Flux (ADU) | Brt px (ADU) | Brt px frac | Top-4 frac |"
        )?;
        writeln!(
            file,
            "|---|---|---|----------:|----------:|-------------:|------------:|-----------:|"
        )?;
        for (i, s) in self.detected_stars.iter().enumerate() {
            writeln!(
                file,
                "| {} | {} | {} | {:.0} | {:.0} | {:.0} | {:.4} | {:.4} |",
                i + 1,
                s.pos.x,
                s.pos.y,
                s.magnitude,
                s.magnitude_adu,
                s.brightest_pixel_adu,
                s.brightest_pixel_part,
                s.top_4_pixels_part
            )?;
        }
        Ok(())
    }

    pub fn width(&self) -> usize {
        self.width as usize
    }

    pub fn height(&self) -> usize {
        self.height as usize
    }

    pub fn quality(&self) -> f64 {
        self.quality
    }

    pub fn exp_t(&self) -> f64 {
        self.exp_t
    }

    pub fn fwhm(&self) -> f64 {
        self.fwhm
    }

    pub fn quality_image(&self) -> Option<f64> {
        self.quality_image
    }

    pub fn quality_star_indices(&self) -> Option<&Vec<usize>> {
        self.quality_star_indices.as_ref()
    }

    pub fn star_count(&self) -> usize {
        self.detected_stars.len()
    }

    pub fn background_raw_adu(&self) -> f64 {
        self.background_level_adu as f64
    }

    pub fn stars(&self) -> &Vec<Star> {
        &self.detected_stars
    }

    /// Calculate the quality score for a selected set without changing the
    /// image's independent quality metric.
    pub fn quality_for_star_indices(
        &self,
        indices: &HashSet<usize>,
        config: &AppConfig,
    ) -> Option<f64> {
        let mut quality_sum = 0.0;
        let mut mag_sum = 0.0;
        for (i, s) in self.detected_stars.iter().enumerate() {
            if indices.contains(&i) && s.magnitude > config.min_photons_quality {
                quality_sum += s.magnitude * s.top_4_pixels_part;
                mag_sum += s.magnitude;
            }
        }
        (mag_sum > 0.0).then_some(quality_sum / mag_sum)
    }

    pub fn stars_with_magnitude_between(
        &self,
        min_magnitude: f64,
        max_magnitude: f64,
    ) -> Vec<usize> {
        let start = self
            .detected_stars
            .partition_point(|s| s.magnitude > max_magnitude);
        let end = self
            .detected_stars
            .partition_point(|s| s.magnitude >= min_magnitude);
        self.detected_stars[start..end]
            .iter()
            .enumerate()
            .map(|(i, _)| start + i)
            .collect()
    }

    pub fn recalculate_quality_for_star_indices(
        &mut self,
        indices: &HashSet<usize>,
        config: &AppConfig,
    ) {
        let mut quality_sum = 0.0;
        let mut mag_sum = 0.0;
        let mut used_indices: Vec<usize> = Vec::new();
        for (i, s) in self.detected_stars.iter().enumerate() {
            if indices.contains(&i) && s.magnitude > config.min_photons_quality {
                quality_sum += s.magnitude * s.top_4_pixels_part;
                mag_sum += s.magnitude;
                used_indices.push(i);
            }
        }
        if mag_sum > 0.0 {
            self.quality_image = Some(self.quality);
            self.quality = quality_sum / mag_sum;
            self.recalculate_fwhm();
        }
        self.quality_star_indices = Some(used_indices);
    }

    fn fmt_impl(&self, f: &mut Formatter<'_>, max_stars: usize) -> std::fmt::Result {
        let n = self.detected_stars.len();
        let avg_top4 = if n > 0 {
            self.detected_stars
                .iter()
                .map(|s| s.top_4_pixels_part)
                .sum::<f64>()
                / n as f64
        } else {
            0.0
        };
        let median_top4 = if n > 0 {
            let mut fwhms: Vec<f64> = self
                .detected_stars
                .iter()
                .map(|s| s.top_4_pixels_part)
                .collect();
            fwhms.sort_by(|a, b| a.total_cmp(b));
            if n.is_multiple_of(2) {
                (fwhms[n / 2 - 1] + fwhms[n / 2]) / 2.0
            } else {
                fwhms[n / 2]
            }
        } else {
            0.0
        };

        writeln!(f, "┌─────────────────────────────────────────┐")?;
        writeln!(
            f,
            "│ Image Metadata Summary [{:>5.2} % ]       │",
            self.quality * 100.0
        )?;
        writeln!(f, "├─────────────────────────────────────────┤")?;
        writeln!(
            f,
            "│  Dimensions   : {:>5} x {:<5}           │",
            self.width, self.height
        )?;
        writeln!(f, "│  Exposure     : {:>8.3} s              │", self.exp_t)?;
        writeln!(f, "│  Gain         : {:>8}                │", self.gain)?;
        writeln!(f, "│  ADU/e⁻       : {:>8.4}                │", self.adu_e)?;
        writeln!(f, "├─────────────────────────────────────────┤")?;
        writeln!(
            f,
            "│  Background   : {:>8} ADU            │",
            self.background_level_adu
        )?;
        writeln!(
            f,
            "│  Sigma (sky)  : {:>8} ADU            │",
            self.sigma_adu
        )?;
        writeln!(f, "├─────────────────────────────────────────┤")?;
        writeln!(f, "│  Stars found  : {:>8}                │", n)?;
        writeln!(
            f,
            "│  Avg top4     : {:>8.2} %              │",
            avg_top4 * 100.0
        )?;
        writeln!(
            f,
            "│  Median top4  : {:>8.2} %              │",
            median_top4 * 100.0
        )?;
        if let Some(ref qi) = self.quality_star_indices {
            writeln!(f, "├─────────────────────────────────────────┤")?;
            writeln!(f, "│  Quality from : {:>4} constellation star │", qi.len())?;
        }
        if let Some(qi) = self.quality_image {
            writeln!(f, "│  Quality (img): {:>8.2} %              │", qi * 100.0)?;
        }
        writeln!(
            f,
            "│  QUALITY      : {:>8.2} %              │",
            self.quality * 100.0
        )?;
        writeln!(
            f,
            "│  FWHM      : {:>8.2} %                 │",
            self.fwhm * 100.0
        )?;
        writeln!(f, "└─────────────────────────────────────────┘")?;

        let star_row = |f: &mut Formatter<'_>, i: usize, s: &Star| -> std::fmt::Result {
            writeln!(
                f,
                "│{:>3} │({:>4},{:>4})   │{:>8.2}    │{:>10.0}  │{:>10.0}  │{:>10.4}  │{:>10.4}  │",
                i + 1,
                s.pos.x,
                s.pos.y,
                s.magnitude,
                s.magnitude_adu,
                s.brightest_pixel_adu,
                s.brightest_pixel_part,
                s.top_4_pixels_part
            )
        };
        let header = || {
            [
                "\n┌────┬──────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐",
                "│ #  │   Position   │  Flux(e⁻)  │  Flux(ADU) │ Brt px ADU │ Brt px frc │ Top4 px fr │",
                "├────┼──────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤",
            ]
        };
        let footer = "└────┴──────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘";

        if let Some(ref qi) = self.quality_star_indices
            && !qi.is_empty()
        {
            writeln!(f, "\n  Top {} Constellation Quality Stars", max_stars)?;
            for line in header() {
                writeln!(f, "{}", line)?;
            }
            let shown = qi.iter().take(max_stars);
            for (i, &idx) in shown.enumerate() {
                star_row(f, i, &self.detected_stars[idx])?;
            }
            writeln!(f, "{}", footer)
        } else {
            let mut by_top4: Vec<&Star> = self.detected_stars.iter().collect();
            by_top4.sort_by(|a, b| b.top_4_pixels_part.total_cmp(&a.top_4_pixels_part));
            writeln!(
                f,
                "\n  Top {} Stars — Best Top-4 Pixel Concentration",
                max_stars
            )?;
            for line in header() {
                writeln!(f, "{}", line)?;
            }
            for (i, s) in by_top4.iter().take(max_stars).enumerate() {
                star_row(f, i, s)?;
            }
            writeln!(f, "{}", footer)
        }
    }
}

//

impl Display for AstroImage {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        self.fmt_impl(f, 10)
    }
}
const DIGIT_BITMAPS: [[u8; 15]; 10] = [
    [1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1], // 0
    [0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1], // 1
    [1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1], // 2
    [1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0], // 3
    [1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1], // 4
    [1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1], // 5
    [0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1], // 6
    [1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1], // 7
    [1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1], // 8
    [1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1], // 9
];

fn draw_digit(
    img: &mut image::RgbImage,
    x0: i32,
    y0: i32,
    scale: i32,
    digit: usize,
    color: image::Rgb<u8>,
) {
    let bitmap = &DIGIT_BITMAPS[digit];
    for row in 0..5i32 {
        for col in 0..3i32 {
            if bitmap[(row * 3 + col) as usize] == 1 {
                for dy in 0..scale {
                    for dx in 0..scale {
                        let px = x0 + col * scale + dx;
                        let py = y0 + row * scale + dy;
                        if px >= 0
                            && py >= 0
                            && (px as u32) < img.width()
                            && (py as u32) < img.height()
                        {
                            img.put_pixel(px as u32, py as u32, color);
                        }
                    }
                }
            }
        }
    }
}

fn draw_number_label(img: &mut image::RgbImage, cx: i32, cy: i32, num: usize) {
    let scale = 3i32;
    let digit_w = 3 * scale;
    let digit_h = 5 * scale;
    let gap = scale;
    let margin = scale;

    let digits: Vec<usize> = num
        .to_string()
        .chars()
        .map(|c| (c as u8 - b'0') as usize)
        .collect();

    let total_w = digits.len() as i32 * digit_w + (digits.len() as i32 - 1) * gap + margin * 2;
    let total_h = digit_h + margin * 2;
    let lx = cx + 7;
    let ly = cy - total_h - 5;

    for dy in 0..total_h {
        for dx in 0..total_w {
            let px = lx + dx;
            let py = ly + dy;
            if px >= 0 && py >= 0 && (px as u32) < img.width() && (py as u32) < img.height() {
                img.put_pixel(px as u32, py as u32, image::Rgb([0u8, 0u8, 0u8]));
            }
        }
    }
    let color = match num {
        x if x < 10 => image::Rgb([200u8, 200u8, 0u8]),
        _ => image::Rgb([50u8, 150u8, 0u8]),
    };

    let mut x_offset = lx + margin;
    for &d in &digits {
        draw_digit(img, x_offset, ly + margin, scale, d, color);
        x_offset += digit_w + gap;
    }
}

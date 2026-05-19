use ndarray::{Array, IxDyn};

pub struct Star {
    pub pos: (usize, usize),
    pub magnitude: f64,
    pub magnitude_adu: f64,
    //pub fwhm: f64,
    pub brightest_pixel_adu: f64,
    pub brightest_pixel_part: f64,
    pub top_4_pixels_part: f64,
    pub ill_defined: bool
}

impl Star {


    pub fn new(pos: (usize, usize), data: &Array<i16, IxDyn>, adu_e: f64, background_adu: u16, psf_size: usize) -> Star {
        let brightest_pixel = (data[[pos.0, pos.1]].wrapping_sub(i16::MIN) as u16 - background_adu) as i32;
        let mut brightest_pixels = [0, 0, 0, 0];
        let mut magnitude = 0.0;
        for i in pos.0-psf_size/2..=pos.0+psf_size/2 {
            for j in pos.1-psf_size/2..=pos.1+psf_size/2 {
                let v = data[[i, j]].wrapping_sub(i16::MIN) as u16;
                let pixel_brightness = v as i32 - background_adu as i32;
                if pixel_brightness > brightest_pixels[0]{
                    brightest_pixels[0] = pixel_brightness;
                    brightest_pixels.sort()
                }
                magnitude += pixel_brightness as f64
            }
        }

        let top_4 = brightest_pixels.into_iter().sum::<i32>();

        Star{
            pos,
            magnitude: magnitude / adu_e,
            magnitude_adu: magnitude,
            brightest_pixel_adu: brightest_pixel as f64,
            brightest_pixel_part: brightest_pixel as f64 / magnitude,
            top_4_pixels_part: top_4 as f64 / magnitude,
            ill_defined: brightest_pixels[3] > brightest_pixel
        }
    }
}
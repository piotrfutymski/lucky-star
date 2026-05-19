use std::collections::HashMap;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref GAIN_ADU: HashMap<u32, f64> = {
        let mut m = HashMap::new();
        m.insert(5200, 27.5);
        m.insert(7000, 40.);
        m.insert(9000, 61.5);
        m.insert(12000, 143.);
        m.insert(15000, 250.);
        m
    };
}

pub const MIN_PHOTONS_TO_DETECT_STAR: i32 = 150;
pub const MIN_CENTRAL_PHOTONS_TO_DETECT_STAR: i32 = 12;
pub const PSF_SIZE: usize = 13;
pub const PSF_SIZE_SQR: usize = PSF_SIZE*PSF_SIZE;
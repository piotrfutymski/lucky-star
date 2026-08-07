use crate::astro_image::AstroImage;
use serde::Deserialize;
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use vector2d::Vector2D;

#[derive(Clone)]
pub struct RegisteredStar {
    pub pos: Vector2D<usize>,
    pub magnitude: f64,
    pub use_in_quality: bool,
    pub median_brightness_adu: Option<f64>,
    pub median_brightest_pixel_part: Option<f64>,
}

// Helper struct that maps the flat JSON fields to RegisteredStar.
#[derive(Deserialize)]
struct RegisteredStarJson {
    x: usize,
    y: usize,
    magnitude: f64,
    use_in_quality: bool,
    #[serde(default)]
    #[serde(alias = "median_brightness")]
    median_brightness_adu: Option<f64>,
    #[serde(default)]
    median_brightest_pixel_part: Option<f64>,
}

impl From<RegisteredStarJson> for RegisteredStar {
    fn from(j: RegisteredStarJson) -> Self {
        RegisteredStar {
            pos: Vector2D::new(j.x, j.y),
            magnitude: j.magnitude,
            use_in_quality: j.use_in_quality,
            median_brightness_adu: j.median_brightness_adu,
            median_brightest_pixel_part: j.median_brightest_pixel_part,
        }
    }
}

impl<'de> Deserialize<'de> for RegisteredStar {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        RegisteredStarJson::deserialize(deserializer).map(RegisteredStar::from)
    }
}

pub fn load_stars_from_json(path: impl AsRef<Path>) -> Result<Vec<RegisteredStar>, Box<dyn Error>> {
    let content = std::fs::read_to_string(path)?;
    let stars: Vec<RegisteredStar> = serde_json::from_str(&content)?;
    if !stars.iter().any(|s| s.use_in_quality) {
        return Err("star pattern must contain at least one star with use_in_quality: true".into());
    }
    Ok(stars)
}

pub struct Constellation {
    pub registered_stars: Vec<RegisteredStar>,
    pub found: bool,
    pub star_mapping: HashMap<usize, usize>,
    pub transform: Option<(Vector2D<f32>, f32)>,
    pub position_tolerance_px: f32,
}

impl Constellation {
    pub fn find_in_image(stars: Vec<RegisteredStar>, image: &AstroImage) -> Constellation {
        Self::find_in_image_with_tolerance(stars, image, 7.0)
    }

    pub fn find_in_image_with_tolerance(
        stars: Vec<RegisteredStar>,
        image: &AstroImage,
        position_tolerance_px: f32,
    ) -> Constellation {
        let mut res = Constellation {
            registered_stars: stars,
            found: false,
            star_mapping: Default::default(),
            transform: None,
            position_tolerance_px,
        };
        let mut mapping = HashMap::new();
        let mut transform = None;
        res.try_to_find_at_image(image, &mut mapping, &mut transform);
        res
    }

    fn try_to_find_at_image(
        &mut self,
        image: &AstroImage,
        test_mapping: &mut HashMap<usize, usize>,
        transform: &mut Option<(Vector2D<f32>, f32)>,
    ) {
        if self.found {
            return;
        }
        let idx = test_mapping.len();
        if idx == self.registered_stars.len() {
            if Self::validate_complete_geometry(
                &self.registered_stars,
                image,
                test_mapping,
                *transform,
                self.position_tolerance_px,
            ) {
                self.found = true;
                self.star_mapping = test_mapping.clone();
                self.transform = *transform;
            }
            return;
        }
        let possible_mappings = Self::find_possible_mappings(
            &self.registered_stars[idx],
            image,
            &self.registered_stars,
            test_mapping,
            *transform,
            self.position_tolerance_px,
        );
        for mapping in possible_mappings {
            // A detected star may represent only one registered star.
            if test_mapping.values().any(|&used| used == mapping) {
                continue;
            }
            test_mapping.insert(idx, mapping);
            if test_mapping.len() == 2 {
                *transform = Some((
                    image.stars()[*test_mapping.get(&0).unwrap()].pos.as_f32s()
                        - self.registered_stars[0].pos.as_f32s(),
                    (image.stars()[*test_mapping.get(&1).unwrap()].pos.as_f32s()
                        - image.stars()[*test_mapping.get(&0).unwrap()].pos.as_f32s())
                    .angle()
                        - (self.registered_stars[1].pos.as_f32s()
                            - self.registered_stars[0].pos.as_f32s())
                        .angle(),
                ))
            }
            self.try_to_find_at_image(image, test_mapping, transform);
            if self.found {
                return;
            }
            if test_mapping.len() == 2 {
                *transform = None;
            }
            test_mapping.remove(&idx);
        }
    }

    fn find_possible_mappings(
        star: &RegisteredStar,
        image: &AstroImage,
        stars: &[RegisteredStar],
        current_mapping: &HashMap<usize, usize>,
        transform: Option<(Vector2D<f32>, f32)>,
        tolerance: f32,
    ) -> Vec<usize> {
        if current_mapping.is_empty() {
            Self::find_in_image_by_geometry(star, image)
        } else if current_mapping.len() == 1 {
            let length_to_match = (star.pos.as_f32s() - stars[0].pos.as_f32s()).length();
            let first_star_pos = image
                .stars()
                .get(*current_mapping.values().next().unwrap())
                .unwrap()
                .pos;
            Self::find_in_image_by_geometry(star, image)
                .into_iter()
                .filter(|m| {
                    Self::is_in_range(
                        length_to_match,
                        first_star_pos,
                        image.stars().get(*m).unwrap().pos,
                        tolerance,
                    )
                })
                .collect()
        } else {
            let transform = transform.unwrap();
            let delta = star.pos.as_f32s() - stars[0].pos.as_f32s();
            let length = delta.length();
            let angle = delta.angle() + transform.1;
            let image_star0 = stars[0].pos.as_f32s() + transform.0;
            let calculated_position =
                image_star0 + Vector2D::new(length * angle.cos(), length * angle.sin());
            Self::find_in_image_by_geometry(star, image)
                .into_iter()
                .filter(|m| {
                    let possible_star_pos = image.stars().get(*m).unwrap().pos.as_f32s();
                    (possible_star_pos - calculated_position).length() <= tolerance
                })
                .collect()
        }
    }

    fn find_in_image_by_geometry(star: &RegisteredStar, image: &AstroImage) -> Vec<usize> {
        let reference_flux = star.magnitude.max(1.0);
        let mut candidates: Vec<usize> = (0..image.stars().len()).collect();
        candidates.sort_by(|&a, &b| {
            let flux_distance = |index: usize| {
                let image_flux = image.stars()[index].magnitude.max(1.0);
                (image_flux / reference_flux).ln().abs()
            };
            flux_distance(a).total_cmp(&flux_distance(b))
        });
        candidates
    }

    fn validate_complete_geometry(
        registered: &[RegisteredStar],
        image: &AstroImage,
        mapping: &HashMap<usize, usize>,
        transform: Option<(Vector2D<f32>, f32)>,
        tolerance: f32,
    ) -> bool {
        let Some((_translation, angle)) = transform else {
            return registered.len() <= 1;
        };
        let (sin, cos) = angle.sin_cos();
        let Some(&mapped_zero) = mapping.get(&0) else {
            return false;
        };
        let image_zero = image.stars()[mapped_zero].pos.as_f32s();
        let reference_zero = registered[0].pos.as_f32s();
        for i in 0..registered.len() {
            let Some(&mapped_i) = mapping.get(&i) else {
                return false;
            };
            let p = registered[i].pos.as_f32s() - reference_zero;
            let predicted = Vector2D::new(
                cos * p.x - sin * p.y + image_zero.x,
                sin * p.x + cos * p.y + image_zero.y,
            );
            let actual = image.stars()[mapped_i].pos.as_f32s();
            if (actual - predicted).length() > tolerance {
                return false;
            }
            for j in 0..i {
                let Some(&mapped_j) = mapping.get(&j) else {
                    return false;
                };
                let expected = (registered[i].pos.as_f32s() - registered[j].pos.as_f32s()).length();
                let observed = (image.stars()[mapped_i].pos.as_f32s()
                    - image.stars()[mapped_j].pos.as_f32s())
                .length();
                if (expected - observed).abs() > tolerance {
                    return false;
                }
            }
        }
        true
    }

    fn is_in_range(
        length_to_match: f32,
        p1: Vector2D<usize>,
        p2: Vector2D<usize>,
        tolerance: f32,
    ) -> bool {
        let length = Vector2D::new(p1.x as f32 - p2.x as f32, p1.y as f32 - p2.y as f32).length();
        (length - length_to_match).abs() <= tolerance
    }
}

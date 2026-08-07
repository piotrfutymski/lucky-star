use std::ops::{Add, Div};

pub fn median<T>(sorted: &[T]) -> Option<T>
where
    T: Copy + Add<Output = T> + Div<Output = T> + From<u8>,
{
    let n = sorted.len();

    if n == 0 {
        return None;
    }

    if n % 2 == 0 {
        Some((sorted[n / 2 - 1] + sorted[n / 2]) / T::from(2))
    } else {
        Some(sorted[n / 2])
    }
}

/// Return <= for PartialOrd types.
///
/// Panics when incomparablee.
fn less_eq<T: PartialOrd>(x: &T, y: &T) -> bool {
    use std::cmp::Ordering::{Less, Equal};

    match x.partial_cmp(y) {
        None => panic!(),
        Some(Less) => true,
        Some(Equal) => true,
        _ => false
    }
}

/// ```
/// use vptree::small_median;
/// assert_eq!(small_median(&[1.0]), 0);
/// assert_eq!(small_median(&[1.0, 2.0]), 0);
/// assert_eq!(small_median(&[2.0, 1.0]), 1);
/// assert_eq!(small_median(&[1.0, 2.0, 3.0]), 1);
/// assert_eq!(small_median(&[1.0, 3.0, 2.0]), 2);
/// assert_eq!(small_median(&[3.0, 1.0, 2.0]), 2);
/// assert_eq!(small_median(&[2.0, 1.0, 3.0]), 0);
/// assert_eq!(small_median(&[3.0, 2.0, 1.0]), 1);
/// assert_eq!(small_median(&[2.0, 3.0, 1.0]), 0);
/// ```
pub fn small_median<T: PartialOrd>(arr: &[T]) -> usize {
    match arr.len() {
        1 => 0,
        2 => if less_eq(&arr[0], &arr[1]) { 0 } else { 1 },
        3 => if less_eq(&arr[1], &arr[2]) {
            if less_eq(&arr[1], &arr[0]) {
                if less_eq(&arr[0], &arr[2]) { 0 } else { 2 }
            } else {
                1
            }
        } else {
            if less_eq(&arr[2], &arr[0]) {
                if less_eq(&arr[0], &arr[1]) { 0 } else { 1 }
            } else {
                2
            }
        },
        _ => unimplemented!()
    }
}


#[cfg(test)]
mod tests
{
    use super::small_median;
    fn test_medians() {
        assert_eq!(small_median(&[1, 2, 3]), 1);
    }
}

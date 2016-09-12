use std::cmp::Ordering;
use rand::{thread_rng, Rng, ThreadRng};
use std::fmt::{Display, Debug};

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


fn small_median_by<T, F>(arr: &[T], f: &F) -> usize
    where F: Fn(&T, &T) -> bool {

    match arr.len() {
        1 => 0,
        2 => if f(&arr[0], &arr[1]) { 0 } else { 1 },
        3 => if f(&arr[1], &arr[2]) {
            if f(&arr[1], &arr[0]) {
                if f(&arr[0], &arr[2]) { 0 } else { 2 }
            } else {
                1
            }
        } else {
            if f(&arr[2], &arr[0]) {
                if f(&arr[0], &arr[1]) { 0 } else { 1 }
            } else {
                2
            }
        },
        _ => unimplemented!()
    }
}

/// Partial sort the elements such that the first k elements are all
/// <= the other elements.
///
/// Takes a function f(x, y) which returns true iff x <= y.
pub fn quick_select_by<T: Display + Debug, F>(arr: &mut [T], k: usize, f: &F) where F: Fn(&T, &T) -> bool {
    let n = arr.len();
    if n <= 1 {
        return;
    }
    if n == 2 {
        // order the elements and return
        if f(&arr[1], &arr[0]) {
            arr.swap(0, 1);
        }
        return;
    }

    // Choose a random pivot (the median among the three random elements)
    let mut rng = thread_rng();
    arr.swap(0, rng.gen_range(0, n));
    arr.swap(1, rng.gen_range(1, n));
    arr.swap(2, rng.gen_range(2, n));
    let mid_idx = small_median_by(&arr[0..3], f);
    arr.swap(0, mid_idx);

    let mut small_it = 1;
    let mut large_it = n - 1;
    loop {
        loop {
            if f(&arr[0], &arr[small_it]) {
                break;
            }
            small_it += 1;
        }
        loop {
            if f(&arr[large_it], &arr[0]) {
                break;
            }
            large_it -= 1;
        }
        if large_it <= small_it {
            break;
        }
        arr.swap(small_it, large_it);
    }

    // Move the median to its correct location.
    arr.swap(0, large_it);

    // Recurse on (at most) one side
    if large_it > k {
        quick_select_by(&mut arr[..large_it], k, f);
    } else if k > large_it {
        quick_select_by(&mut arr[(large_it + 1)..], k - large_it - 1, f);
    }
}


#[cfg(test)]
mod tests
{
    use super::small_median;
    use super::quick_select_by;

    #[test]
    fn test_medians() {
        assert_eq!(small_median(&[1, 2, 3]), 1);
    }

    #[test]
    fn test_quick_select_by() {
        for i in 0..11 {
            let mut v = vec![2, 0, 4, 6, 5, 1, 3, 9, 7, 8, 2];
            quick_select_by(&mut v, i, &|a, b| a <= b);
            for x in 0..i {
                assert!(v[x] <= v[i]);
            }
            for x in i+1..10 {
                assert!(v[i] <= v[x]);
            }
        }
    }
}

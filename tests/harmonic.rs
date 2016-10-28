extern crate vptree;

use vptree::{MetricItem, VPTree};

#[derive(Debug, PartialEq, Clone)]
struct Point(f32);

impl MetricItem<f32> for Point {
    fn distance(&self, a: &Self) -> f32 {
        return (self.0 - a.0).abs()
    }
}


#[test]
fn test_linear() {
    let mut done = false;
    for n in 10..101 {
        let points: Vec<_> = (0..n+1).map(|x| Point(x as f32/ n as f32)).collect();

        let vp = VPTree::new(points.iter().cloned().collect());
        assert!(vp.is_some());
        let vp = vp.unwrap();

        for p in &points {
            // each point should be nearest to itself.
            if vp.nearest_neighbor(p) != p {
                print!("\n{}", vp.dump());
                done = true;
            }
            assert_eq!(vp.nearest_neighbor(p), p);
         }
        if done {
            break;
        }
    }
}

#[test]
fn test_harmonic() {
    let mut done = false;
    for n in 10..101 {
        let points: Vec<_> = (1..n+1).map(|x| Point(1.0 / (x as f32))).collect();

        let vp = VPTree::new(points.iter().cloned().collect());
        assert!(vp.is_some());
        let vp = vp.unwrap();

        for p in &points {
            // each point should be nearest to itself.
            if vp.nearest_neighbor(p) != p {
                print!("\n{}", vp.dump());
                done = true;
            }
            assert_eq!(vp.nearest_neighbor(p), p);
         }
        if done {
            break;
        }
    }
}

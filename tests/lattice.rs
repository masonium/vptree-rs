extern crate vptree;

use vptree::{MetricItem, VPTree};

#[derive(Debug)]
struct Point {
    x: f32,
    y: f32
}
impl Point {
    fn new(x: f32, y: f32) -> Self {
        Point { x: x, y: y }
    }
}

impl MetricItem<f32> for Point {
    fn distance(&self, q: &Self) -> f32 {
        let dx = self.x - q.x;
        let dy = self.y - q.y;
        (dx*dx + dy*dy).sqrt()
    }
}

#[test]
fn point_check() {
    let a = Point::new(0.0, 0.0);
    let b = Point::new(1.0, 0.0);

    assert_eq!(Point::distance(&a, &b), 1.0);
}

fn lattice_points(n: usize) -> Vec<Point> {
    (0..n).flat_map( |i| {
        (0..n).map(move |j| {
            Point::new(i as f32, j as f32)
        })
    }).collect()
}

#[test]
fn lattice_vpn() {
    let points: Vec<Point> = lattice_points(20);

    let tree = VPTree::new(points).unwrap();

    let ps = tree.nearest_neighbors(&Point::new(4.46, 4.4), 4, true);
    assert_eq!(ps.len(), 4);
    assert_eq!(ps[0].x, 4.0);
    assert_eq!(ps[0].y, 4.0);

    assert_eq!(ps[1].x, 5.0);
    assert_eq!(ps[1].y, 4.0);

    assert_eq!(ps[2].x, 4.0);
    assert_eq!(ps[2].y, 5.0);

    assert_eq!(ps[3].x, 5.0);
    assert_eq!(ps[3].y, 5.0);
}

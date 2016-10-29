//! 'vptree-rs` is a crate for vantage point trees.
//!
//! A vantage point tree, or VP-tree is a metric tree data structure for
//! nearest neighbor and range queries. A VP-tree stores a set of points
//! related by a _metric_, a function f with the following properties:
//!
//! For any two distinct points x and y:
//!
//! - f(x, x) = 0
//! - f(x, y) > 0
//! - f(x, y) = f(y, x)
//! - f(x, z) <= f(x, y) + f(y, z)
//!
//! `VPTree`s in the `vptree-rs` crate support k-nearest neighbor and
//! radius queries.
//!
//! VP-trees work by recursively splitting the data set in two, based on
//! how far each point is from a selected _vantage point_. When searching
//! through the VP-tree, we use the distance from query point to a
//! subtree's vantage point to potentially cull a subtree.
//!
//! # Examples
//!
//! ```rust
//! use self::vptree::{MetricItem, VPTree};
//!
//! #[derive(Debug)]
//! struct Point {
//!     x: f32,
//!     y: f32
//! }
//! impl Point {
//!     fn new(x: f32, y: f32) -> Self {
//!         Point { x: x, y: y }
//!     }
//! }
//!
//! impl MetricItem<f32> for Point {
//!     fn distance(&self, q: &Self) -> f32 {
//!         let dx = self.x - q.x;
//!         let dy = self.y - q.y;
//!         (dx*dx + dy*dy).sqrt()
//!     }
//! }
//!
//! fn lattice_points(n: usize) -> Vec<Point> {
//!     (0..n).flat_map( |i| {
//!         (0..n).map(move |j| {
//!             Point::new(i as f32, j as f32)
//!         })
//!     }).collect()
//! }
//!
//! #[test]
//! fn lattice_vpn() {
//!     let points: Vec<Point> = lattice_points(10);
//!
//!     let tree = VPTree::new(points).unwrap();
//!
//!     let ps = tree.nearest_neighbors(&Point::new(4.46, 4.4), 4, true);
//!     assert_eq!(ps.len(), 4);
//!     assert_eq!(ps[0].x, 4.0);
//!     assert_eq!(ps[0].y, 4.0);
//!
//!     assert_eq!(ps[1].x, 5.0);
//!     assert_eq!(ps[1].y, 4.0);
//!
//!     assert_eq!(ps[2].x, 4.0);
//!     assert_eq!(ps[2].y, 5.0);
//!
//!     assert_eq!(ps[3].x, 5.0);
//!     assert_eq!(ps[3].y, 5.0);
//! }
//! ```
//!

extern crate num;
extern crate rand;
extern crate order_stat;

pub mod vptree;

pub use vptree::{VPTree, MetricItem};

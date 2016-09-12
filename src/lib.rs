extern crate num;
extern crate rand;

pub mod vptree;
pub mod median;

pub use vptree::{VPTree, MetricItem};
pub use median::small_median;

//! Vantage-Point Trees are a data structure for fast
//! k-nearest-neighbor searches.
extern crate rand;

use rand::distributions::{Range, IndependentSample};
use std::borrow::Borrow;
use std::collections::{BinaryHeap};
use std::cmp::{Ord, PartialOrd, Ordering};
use std::fmt::{Debug, Display};
pub use num::Float;

pub trait BoundedMetricItem<F: Float> {
    //// Return the standard bounded-distance transformation, for a
    //// given distance.
    fn bounded_distance(a: &Self, b: &Self) -> F;
}

// impl MetricItem for f32 {
//     type Measure = f32;
// }

struct TaggedItem<F: Float, T: BoundedMetricItem<F>> {
    pub item: T,
    hist: Vec<F>
}

/// Randomly select a vantage point, and choose the point furthers
/// from that point.
fn select_vantage_point<F: Float, T: BoundedMetricItem<F>>(items: &Vec<TaggedItem<F, T>>) -> usize {
    let mut rng = rand::thread_rng();

    let range = Range::new(0, items.len());
    let i = range.ind_sample(&mut rng);
    let random_item = &items[i];

    let min_d = (T::bounded_distance(&random_item.item, &random_item.item), i);

    items.iter().enumerate().fold(min_d, |acc, (i, y)| {
        let d = T::bounded_distance(&random_item.item, &y.item);
        if d < acc.0 { (d, i) } else { acc }
    }).1
}

pub trait Scalar  : Float + Debug + Display {}
impl<T: Float + Debug + Display> Scalar for T {}

struct VPNode<F: Scalar, T: BoundedMetricItem<F>> {
    inner: Option<Box<VPNode<F, T>>>,
    outer: Option<Box<VPNode<F, T>>>,
    center: T,
    mu: Option<F>
}

struct HeapElem<'a, F: Scalar, T: 'a> {
    dist: F,
    item: &'a T
}

impl<'a, F: Scalar, T: 'a> HeapElem<'a, F, T> {
    fn new(d: F, i: &'a T) -> Self{
        HeapElem { dist: d, item: i }
    }
}

impl<'a, F: Scalar, T: 'a> PartialOrd for HeapElem<'a, F, T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.dist.partial_cmp(&other.dist)
    }
}


impl<'a, F: Scalar, T: 'a> PartialEq for HeapElem<'a, F, T> {
    fn eq(&self, other: &Self) -> bool {
        self.dist.eq(&other.dist)
    }
}

impl<'a, F: Scalar, T: 'a> Eq for HeapElem<'a, F, T> {
}

impl<'a, F: Scalar, T: 'a> Ord for HeapElem<'a, F, T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

// impl<'a, F: Scalar, T: 'a> PartialOrd for HeapElem<'a, F, T> {
//     fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
//         self.dist.parital_cmp(other.dist);
//     }
// }


impl<F: Scalar, T: BoundedMetricItem<F>> VPNode<F, T> {
    // new
    pub fn new(mut items: Vec<TaggedItem<F, T>>) -> VPNode<F, T> {
        if items.len() == 1 {
            return VPNode { inner: None,
                            outer: None,
                            center: items.pop().unwrap().item,
                            mu: None };
        }

        let sel_index = select_vantage_point(&items);

        let n = items.len();
        let vp = items.swap_remove(sel_index);

        // sort the items by distance from the selected item.
        for mut ti in items.iter_mut() {
            let d = T::bounded_distance(&ti.item, &vp.item);
            ti.hist.push(d);
        }

        // split the items into sectinos, where one section is less than the other
        items.sort_by(|a, b| a.hist.last().unwrap().partial_cmp(b.hist.last().unwrap()).unwrap());

        let right_items = items.split_off((n + 1) / 2);
        let mu = *(items.last().unwrap().hist.last().unwrap());
        let inner = if items.is_empty() { None } else { Some(Box::new(VPNode::new(items))) };
        let outer = if right_items.is_empty() { None } else { Some(Box::new(VPNode::new(right_items))) };

        VPNode { inner: inner, outer: outer, center: vp.item, mu: Some(mu)  }
    }

    /// Push the nearest neighbors of this tree onto the binary heap,
    /// replacing existing further-away elemtns as necessary.
    pub fn nearest_neighbors<'b, 'a: 'b>(&'a self, obj: &T, n: usize, traversed: &mut usize,
                                         heap: &'b mut BinaryHeap<HeapElem<'a, F, Self>>)  {
        *traversed += 1;
        let d_center = T::bounded_distance(obj, &self.center);

        let elem = HeapElem::new(d_center, self);

        // Push the element on if it is closer than the current furthest element.
        if heap.len() < n {
            heap.push(elem);
        } else if heap.peek().unwrap().dist > elem.dist {
            heap.pop();
            heap.push(elem);
        }

        // If we have an inner or outer node.
        if let Some(mu) = self.mu {
            let mut nodes = [(&self.inner, true), (&self.outer, false)];

            // Traverse the outer node first if we're outside the ring.
            if d_center > mu {
                nodes.swap(0, 1);
            }

            for &(node_opt, is_inner) in &nodes {
                if let Some(ref node) = *node_opt {
                    let d_max = heap.peek().unwrap().dist;
                    let possible_new_elem = (is_inner && d_max > d_center - mu) || (!is_inner && d_max > mu - d_center);
                    if possible_new_elem {
                        let x: &Self = node.borrow();
                        x.nearest_neighbors(obj, n, traversed, heap);
                    }
                }
            }
        }
    }
}

/// implement a Vp-s tree
pub struct VPTree<F: Scalar, T: BoundedMetricItem<F>> {
    root: VPNode<F, T>
}

impl<F: Scalar, T: BoundedMetricItem<F>> VPTree<F, T> {
    /// Construct a new vantage point tree from a set of elements.
    pub fn new(items: Vec<T>) -> Option<VPTree<F, T>> {
        let n = items.len();
        let ln = (n as f64).log2().ceil() as usize;
        if n > 0 {
            let tagged_items: Vec<TaggedItem<F, T>> = items.into_iter()
                .map(|x| TaggedItem { item: x, hist: Vec::with_capacity(ln) }).collect();
            Some(VPTree { root: VPNode::new(tagged_items) })
        } else {
            None
        }
    }

    /// find the nearest neighbor
    pub fn nearest_neighbor(&self, obj: &T) -> &T {
        let mut heap = BinaryHeap::with_capacity(1);
        let mut traversed: usize = 0;
        self.root.nearest_neighbors(obj, 1, &mut traversed, &mut heap);

        let he = heap.pop().unwrap();
        &he.item.center
    }

    /// Find the n nearest neighbors.
    pub fn nearest_neighbors(&self, obj: &T, n: usize, traversed: &mut usize) -> Vec<&T> {
        let mut heap = BinaryHeap::with_capacity(n);
        *traversed = 0;
        self.root.nearest_neighbors(obj, n, traversed, &mut heap);

        heap.into_sorted_vec().into_iter().map(|x| &x.item.center).collect()
    }
}

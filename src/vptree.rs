//! Vantage-Point Trees are a data structure for fast
//! k-nearest-neighbor searches.
extern crate rand;

use rand::distributions::{Range, IndependentSample};
use std::borrow::Borrow;
use std::collections::{BinaryHeap};
use std::cmp::{Ord, PartialOrd, Ordering};
pub use num::Float;


pub trait MetricItem<F: Float> {
    fn distance(a: &Self, b: &Self) -> F;
}

pub trait BoundedMetricItem<F: Float> {
    //// Return the standard bounded-distance transformation, for a
    //// given distance.
    fn bounded_distance(a: &Self, b: &Self) -> F;
}

impl<X: MetricItem<f32>> BoundedMetricItem<f32> for X  {
    fn bounded_distance(a: &Self, b: &Self) -> f32 {
        let d = X::distance(a, b);
        d / (1.0 + d)
    }
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

struct VPNode<F: Float, T: BoundedMetricItem<F>> {
    inner: Option<Box<VPNode<F, T>>>,
    outer: Option<Box<VPNode<F, T>>>,
    center: T,
    mu: Option<F>
}

struct HeapElem<'a, F: Float, T: 'a> {
    dist: F,
    item: &'a T
}

impl<'a, F: Float, T: 'a> HeapElem<'a, F, T> {
    fn new(d: F, i: &'a T) -> Self{
        HeapElem { dist: d, item: i }
    }
}

impl<'a, F: Float, T: 'a> PartialOrd for HeapElem<'a, F, T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.dist.partial_cmp(&other.dist)
    }
}


impl<'a, F: Float, T: 'a> PartialEq for HeapElem<'a, F, T> {
    fn eq(&self, other: &Self) -> bool {
        self.dist.eq(&other.dist)
    }
}

impl<'a, F: Float, T: 'a> Eq for HeapElem<'a, F, T> {
}

impl<'a, F: Float, T: 'a> Ord for HeapElem<'a, F, T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

// impl<'a, F: Float, T: 'a> PartialOrd for HeapElem<'a, F, T> {
//     fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
//         self.dist.parital_cmp(other.dist);
//     }
// }


/// Push elements on to a max heap, preserving the n-closest
/// elements.
fn push_up_to_max<'a, 'b: 'a, F: Float, T: BoundedMetricItem<F>>(heap: &'b mut BinaryHeap<HeapElem<'a, F, T>>,
                                                                 n: usize,
                                                                 elem: HeapElem<'a, F, T>) {
    if heap.len() < n {
        heap.push(elem);
    } else if heap.peek().unwrap().dist > elem.dist {
        heap.replace(elem);
    }
}

impl<F: Float, T: BoundedMetricItem<F>> VPNode<F, T> {
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

        // split the items into
        items.sort_by(|a, b| a.hist.last().unwrap().partial_cmp(b.hist.last().unwrap()).unwrap());

        let right_items = items.split_off((n + 1) / 2);
        let mu = *(items.last().unwrap().hist.last().unwrap());
        let inner = if items.is_empty() { None } else { Some(Box::new(VPNode::new(items))) };
        let outer = if right_items.is_empty() { None } else { Some(Box::new(VPNode::new(right_items))) };

        VPNode { inner: inner, outer: outer, center: vp.item, mu: Some(mu)  }
    }

    pub fn nearest_neighbors<'b, 'a: 'b>(&'a self, obj: &T, n: usize, heap: &'b mut BinaryHeap<HeapElem<'a, F, Self>>) {
        let d = T::bounded_distance(obj, &self.center);

        let elem = HeapElem::new(d, self);
        if heap.len() < n {
            heap.push(elem);
        } else if heap.peek().unwrap().dist > elem.dist {
            heap.replace(elem);
        }

        if let Some(ref inner) = self.inner {
            let x: &Self = inner.borrow();
            x.nearest_neighbors(obj, n, heap);
        }
        if let Some(ref outer) = self.outer {
            let x: &Self = outer.borrow();
            x.nearest_neighbors(obj, n, heap);
        }
    }

    /// Return the closest element to the input element.
    pub fn nearest_neighbor(&self, obj: &T) -> (&T, F) {
        unimplemented!();
        // let d = T::bounded_distance(obj, &self.center);
        // let pair = (&self.center, d);
        // match self.mu {
        //     None => pair,
        //     Some(dist) => {
        //         // If we're closer to the center than the edge, we
        //         // know we don't have to look on the outer section.
        //         if d * (F::one() + F::one())  < dist {
        //             match self.inner {
        //                 Some(ref inner) => {
        //                     let node: &VPNode<_, _> = inner.borrow();
        //                     let inner_nn = node.nearest_neighbor(obj);
        //                     if inner_nn.1 < d { inner_nn } else { pair }
        //                 },
        //                 None => pair
        //             }
        //         } else {
        //             pair
        //         }
        //     }
        // }
    }
}

/// implement a Vp-s tree
pub struct VPTree<F: Float, T: BoundedMetricItem<F>> {
    root: VPNode<F, T>
}

impl<F: Float, T: BoundedMetricItem<F>> VPTree<F, T> {
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
    pub fn nearest_neighbor(&self, obj: &T) -> (&T, F) {
        self.root.nearest_neighbor(obj)
    }
}

/// Re-arrange the list such that item k is k-th in the list, and all items > k are
fn quick_select<T: PartialOrd + Copy>(s: &mut [T], k: usize) -> usize {
    // If the range is smaller than some limit, brute-force the median.
    unimplemented!();
}

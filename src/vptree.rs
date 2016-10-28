//! Vantage-Point Trees are a data structure for fast
//! k-nearest-neighbor searches.
extern crate rand;

use rand::distributions::{Range, IndependentSample};
use std::borrow::Borrow;
use std::collections::{BinaryHeap};
use std::cmp::{Ord, PartialOrd, Ordering};
use std::fmt::{Debug, Display};
use std::ops::Mul;
pub use num::Float;
use order_stat::kth_by;

pub trait MetricValue: Float + Mul<f32, Output=Self> {
}
impl<T: Float + Mul<f32, Output=T>> MetricValue for T {
}

/// Defines a metric for items in a metric space.
///
/// A metric is a function on a set S, with the following properties.
/// - f(x, y) >= 0, with equality iff x == y
/// - f(x, y) = f(y, x)
/// - f(x, y) + f(y, z) <= f(x, z)
///
/// A VP-Tree can only be constructed by a set forming a metric. If
/// the `distance` function does not satisfy the metric conditions, a
/// vp-tree constructed from the elements will not be correct.
pub trait MetricItem<F: MetricValue> {
    /// Return the distance to another element in the metric space.
    ///
    /// The `distance` function must satisfy the metric properties.
    /// (See: )
    fn distance(&self, b: &Self) -> F;
}

struct TaggedItem<F: MetricValue, T: MetricItem<F>> {
    pub item: T,
    pub dist: F
}

/// Return a randomly-selected vantage point.
fn select_vantage_point<F: MetricValue, T: MetricItem<F>>(items: &Vec<TaggedItem<F, T>>) -> usize {
    // Randomly select a point.
    let mut rng = rand::thread_rng();

    let range = Range::new(0, items.len());
    let i = range.ind_sample(&mut rng);
    let random_item = &items[i];

    let min_d = (F::zero(), i);

    // The vantage point will be the point furthest from the selected
    // one.
    items.iter().enumerate().fold(min_d, |acc, (i, y)| {
        let d = T::distance(&random_item.item, &y.item);
        if d > acc.0 { (d, i) } else { acc }
    }).1
}

pub trait Scalar : MetricValue + Debug + Display {}
impl<T: MetricValue + Debug + Display> Scalar for T {}

struct InnerNode<F: Float, N> {
    pub mu: F,
    pub inner: Box<N>,
    pub outer: Option<Box<N>>
}

struct VPNode<F: Scalar, T: MetricItem<F>> {
    contents: Option<InnerNode<F, VPNode<F, T>>>,
    center: T,
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

impl<F: Scalar, T: MetricItem<F>> VPNode<F, T> {
    // new
    pub fn new(mut items: Vec<TaggedItem<F, T>>) -> VPNode<F, T> {
        if items.len() == 1 {
            return VPNode { contents: None,
                            center: items.pop().unwrap().item
            };
        }

        let sel_index = select_vantage_point(&items);

        let vp = items.swap_remove(sel_index);

        // Compute the new distance from the vantage point for all of
        // the items.
        for mut ti in items.iter_mut() {
            ti.dist = T::distance(&ti.item, &vp.item);
        }

        let n = items.len();

        // We want to split the array into two as follows:
        //
        // The left array gets an extra element when the number of
        // elements is odd.
        //
        // The last element of the left array is larger than all
        // others, and smaller than eevery element in the right array.
        if n > 1 {
            kth_by(&mut items, (n-1)/2, |a, b| a.dist.partial_cmp(&b.dist).unwrap());
        }

        let right_items = items.split_off((n+1)/2);

        match items.last().map(|x| x.dist) {
            Some(dist) => {
                let mu = dist;
                let inner = Box::new(VPNode::new(items));
                let outer = if right_items.is_empty() { None } else { Some(Box::new(VPNode::new(right_items))) };
                VPNode { center: vp.item, contents: Some(InnerNode {
                    mu: mu,
                    inner: inner,
                    outer: outer })}
            },
            None => {
                VPNode { center: vp.item, contents: None }
            }
        }
    }

    /// Push the nearest neighbors of this tree onto the binary heap,
    /// replacing existing further-away elemtns as necessary.
    pub fn nearest_neighbors<'b, 'a: 'b>(&'a self, obj: &T, n: usize,
                                         heap: &'b mut BinaryHeap<HeapElem<'a, F, Self>>)  {
        let d_center = T::distance(obj, &self.center);

        let elem = HeapElem::new(d_center, self);

        // Push the element on if it is closer than the current furthest element.
        if heap.len() < n {
            heap.push(elem);
        } else if heap.peek().unwrap().dist > elem.dist {
            heap.pop();
            heap.push(elem);
        }

        // If we have an inner or outer node.
        if let Some(ref contents) = self.contents {
            let mu = contents.mu;
            let some_inner = Some(&contents.inner);
            let mut nodes = [(some_inner, true), (contents.outer.as_ref(), false)];

            // Traverse the outer node first if we're outside the ring.
            if d_center > mu {
                nodes.swap(0, 1);
            }

            for &(node_opt, is_inner) in &nodes {
                if let Some(node) = node_opt {
                    let d_max = heap.peek().unwrap().dist;
                    let possible_new_elem = (is_inner && d_max > d_center - mu) || (!is_inner && d_max > mu - d_center);
                    if possible_new_elem {
                        let x: &Self = node.borrow();
                        x.nearest_neighbors(obj, n, heap);
                    }
                }
            }
        }
    }

    /// Return all elements within a given radius of the node.
    pub fn within_radius<'a, 'b: 'a>(&'b self, obj: &T, radius: F, v: &mut Vec<HeapElem<'a, F, Self>>) {
        let d_center = T::distance(obj, &self.center);

        // Push the element on if it is closer than the current furthest element.
        if d_center < radius {
            v.push(HeapElem::new(d_center, &self));
        }

        // If we have an inner or outer node.
        if let Some(ref contents) = self.contents {
            let mu = contents.mu;
            let some_inner = Some(&contents.inner);
            let mut nodes = [(some_inner, true), (contents.outer.as_ref(), false)];

            // Traverse the outer node first if we're outside the ring.
            if d_center > mu {
                nodes.swap(0, 1);
            }

            for &(node_opt, is_inner) in &nodes {
                if let Some(node) = node_opt {
                    let possible_new_elem = (is_inner && radius > d_center - mu) || (!is_inner && radius > mu - d_center);
                    if possible_new_elem {
                        let x: &Self = node.borrow();
                        x.within_radius(obj, radius, v);
                    }
                }
            }
        }

    }
}

/// implement a Vp-s tree
pub struct VPTree<F: Scalar, T: MetricItem<F>> {
    root: VPNode<F, T>
}

impl<F: Scalar, T: MetricItem<F>> VPTree<F, T> {
    /// Construct a new vantage point tree from a set of elements.
    pub fn new(items: Vec<T>) -> Option<VPTree<F, T>> {
        let n = items.len();
        if n > 0 {
            let tagged_items: Vec<TaggedItem<F, T>> = items.into_iter()
                .map(|x| TaggedItem { item: x, dist: F::zero() }).collect();
            Some(VPTree { root: VPNode::new(tagged_items) })
        } else {
            None
        }
    }

    /// Return all elements with a given radius of the target.
    ///
    /// If `sorted` is true, the elements are sorted by ascending
    /// distance from the query point,
    pub fn within_radius(&self, obj: &T, radius: F, sorted: bool) -> Vec<&T> {
        let mut elems = Vec::new();
        self.root.within_radius(obj, radius, &mut elems);

        if sorted {
            elems.sort();
        }

        elems.into_iter().map(|x| &x.item.center).collect()
    }

    /// Find the nearest neighbor.
    pub fn nearest_neighbor(&self, obj: &T) -> &T {
        let mut heap = BinaryHeap::with_capacity(1);
        self.root.nearest_neighbors(obj, 1, &mut heap);

        let he = heap.pop().unwrap();
        &he.item.center
    }

    /// Find the n nearest neighbors. IF sorted is true, the results
    /// will be sorted from nearest to furthest.
    pub fn nearest_neighbors(&self, obj: &T, n: usize, sorted: bool) -> Vec<&T> {
        let mut heap = BinaryHeap::with_capacity(n);
        self.root.nearest_neighbors(obj, n, &mut heap);

        let v = if sorted {
            heap.into_sorted_vec()
        } else {
            heap.into_vec()
        };
        v.into_iter().map(|x| &x.item.center).collect()

    }
}
impl<F: Scalar, T: MetricItem<F> + Debug> VPNode<F, T> {
    pub fn dump(&self, prefix: &str) -> String {
        let mut s: String = format!("{}elem: {:?}", prefix, self.center);
        if let Some(ref c) = self.contents {
            s += &format!(", mu: {}\n", c.mu);
            let new_prefix = format!("{}  ", prefix);

            let ref n: VPNode<F, T> = *c.inner.borrow();
            s += &format!("{}{}", prefix, n.dump(&new_prefix));

            if let Some(ref outer) = c.outer {
                let ref n: VPNode<F, T> = *outer.borrow();
                s += &format!("{}{}", prefix, n.dump(&new_prefix));
            }
        }

        s

    }
}

impl <F: Scalar, T: MetricItem<F> + Debug> VPTree<F, T> {
    pub fn dump(&self) -> String {
        self.root.dump("")
    }
}

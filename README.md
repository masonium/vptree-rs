# vptree-rs [![Build Status](https://travis-ci.org/masonium/vptree-rs.svg?branch=master)](https://travis-ci.org/masonium/vptree-rs) #
`vptree-rs`is an implementation of vantage point trees in rust.

Vantage Point trees are a metric tree data structure, used for
nearest-neighbor or range queries. It is constructed from a fixed set
of points. You can then, for any novel point, find the exact k-nearest
neighbors or all points within a given radius.

## Example
TODO: fill in example

## Limitations
- Does not yet support a combined `nearest_neighbors_within_range`
  query.

//! Copied over from a #654 in itertools until that is merged.
//! https://github.com/rust-itertools/itertools/pull/654
use core::cmp::{Ordering, Reverse};
use std::vec::{IntoIter, Vec};

/// Consumes a given iterator, returning the minimum elements in **ascending** order.
fn k_smallest_general<T, I: Iterator<Item = T>>(
    mut iter: I,
    k: usize,
    mut comparator: impl FnMut(&T, &T) -> Ordering,
) -> IntoIter<T> {
    if k == 0 {
        return Vec::new().into_iter();
    }
    let mut storage: Vec<T> = iter.by_ref().take(k).collect();
    let mut heap = &mut storage[..];

    // Rearrange the into a valid heap by reordering from the second-bottom-most layer up
    // Slightly faster than ordering on each insert
    // (But only by a factor of lg(k) and I'd love to hear of a use case where that matters)
    for i in (0..(heap.len() / 2 + 1)).rev() {
        sift_down(heap, &mut comparator, i);
    }

    if k == heap.len() {
        // Nothing else needs done if we didn't fill the storage in the first place
        // Also avoids unexpected behaviour with restartable iterators
        for val in iter {
            if comparator(&heap[0], &val) == Ordering::Greater {
                heap[0] = val;
                sift_down(heap, &mut comparator, 0);
            }
        }
    }

    while heap.len() > 1 {
        let last_idx = heap.len() - 1;
        heap.swap(0, last_idx);
        // Leaves the length shorter than the number of elements
        // so that sifting does not disturb already popped elements
        heap = &mut heap[..last_idx];
        sift_down(heap, &mut comparator, 0);
    }

    storage.into_iter()
}

fn reverse_cmp<T, F>(cmp: F) -> impl Fn(&T, &T) -> Ordering
where
    F: Fn(&T, &T) -> Ordering,
{
    move |a, b| cmp(b, a)
}

fn reverse_key<T, K, F>(f: F) -> impl Fn(&T) -> Reverse<K>
where
    F: Fn(&T) -> K,
{
    move |a| Reverse(f(a))
}

/// Sift the element currently at `origin` **away** from the root until it is properly ordered
fn sift_down<T>(
    heap: &mut [T],
    comparator: &mut impl FnMut(&T, &T) -> Ordering,
    mut origin: usize,
) {
    fn children_of(n: usize) -> (usize, usize) {
        (2 * n + 1, 2 * n + 2)
    }

    while origin < heap.len() {
        let (left_idx, right_idx) = children_of(origin);
        if left_idx >= heap.len() {
            return;
        }

        let replacement_idx = if right_idx < heap.len()
            && Ordering::Less == comparator(&heap[left_idx], &heap[right_idx])
        {
            right_idx
        } else {
            left_idx
        };

        let cmp = comparator(&heap[origin], &heap[replacement_idx]);
        if Ordering::Less == cmp {
            heap.swap(origin, replacement_idx);
            origin = replacement_idx;
        } else {
            return;
        }
    }
}

pub trait ItertoolsExt: Iterator {
    /// itertools::assert_equal(five_smallest, 0..5);
    /// ```
    fn k_smallest(mut self, k: usize) -> std::vec::IntoIter<Self::Item>
    where
        Self: Sized,
        Self::Item: Ord,
    {
        // The stdlib heap has optimised handling of "holes", which is not included in our heap implementation in k_smallest_general.
        // While the difference is unlikely to have practical impact unless `T` is very large, this method uses the stdlib structure
        // to maintain performance compared to previous versions of the crate.
        use std::collections::BinaryHeap;

        if k == 0 {
            return vec![].into_iter();
        }

        let mut heap = self.by_ref().take(k).collect::<BinaryHeap<_>>();

        self.for_each(|i| {
            debug_assert_eq!(heap.len(), k);
            // Equivalent to heap.push(min(i, heap.pop())) but more efficient.
            // This should be done with a single `.peek_mut().unwrap()` but
            //  `PeekMut` sifts-down unconditionally on Rust 1.46.0 and prior.
            if *heap.peek().unwrap() > i {
                *heap.peek_mut().unwrap() = i;
            }
        });

        heap.into_sorted_vec().into_iter()
    }

    /// Sort the k smallest elements into a new iterator using the provided comparison.
    ///
    /// This corresponds to `self.sorted_by(cmp).take(k)` in the same way that
    /// [Itertools::k_smallest] corresponds to `self.sorted().take(k)`, in both semantics and complexity.
    /// Particularly, a custom heap implementation ensures the comparison is not cloned.
    fn k_smallest_by<F>(self, k: usize, cmp: F) -> std::vec::IntoIter<Self::Item>
    where
        Self: Sized,
        F: Fn(&Self::Item, &Self::Item) -> Ordering,
    {
        k_smallest_general(self, k, cmp)
    }

    /// Return the elements producing the k smallest outputs of the provided function
    ///
    /// This corresponds to `self.sorted_by_key(cmp).take(k)` in the same way that
    /// [Itertools::k_smallest] corresponds to `self.sorted().take(k)`, in both semantics and time complexity.
    /// This method will use an _additional_ `k * sizeof(K)` memory compared to that method.
    fn k_smallest_by_key<F, K>(self, k: usize, key: F) -> std::vec::IntoIter<Self::Item>
    where
        Self: Sized,
        F: Fn(&Self::Item) -> K,
        K: Ord,
    {
        self.k_smallest_by(k, |a, b| key(a).cmp(&key(b)))
    }

    /// Sort the k largest elements into a new iterator, in descending order.
    /// Semantically equivalent to `k_smallest` with a reversed `Ord`
    /// However, this is implemented by way of a custom binary heap
    /// which does not have the same performance characteristics for very large `T`
    /// ```
    /// use itertools::Itertools;
    ///
    /// // A random permutation of 0..15
    /// let numbers = vec![6, 9, 1, 14, 0, 4, 8, 7, 11, 2, 10, 3, 13, 12, 5];
    ///
    /// let five_largest = numbers
    ///     .into_iter()
    ///     .k_largest(5);
    ///
    /// itertools::assert_equal(five_largest, vec![14,13,12,11,10]);
    /// ```
    fn k_largest(self, k: usize) -> std::vec::IntoIter<Self::Item>
    where
        Self: Sized,
        Self::Item: Ord,
    {
        self.k_smallest_by(k, reverse_cmp(Self::Item::cmp))
    }

    /// Sort the k largest elements into a new iterator using the provided comparison.
    /// Functionally equivalent to `k_smallest_by` with a reversed `Ord`
    fn k_largest_by<F>(self, k: usize, cmp: F) -> std::vec::IntoIter<Self::Item>
    where
        Self: Sized,
        F: Fn(&Self::Item, &Self::Item) -> Ordering,
    {
        self.k_smallest_by(k, reverse_cmp(cmp))
    }

    /// Return the elements producing the k largest outputs of the provided function
    fn k_largest_by_key<F, K>(self, k: usize, key: F) -> std::vec::IntoIter<Self::Item>
    where
        Self: Sized,
        F: Fn(&Self::Item) -> K,
        K: Ord,
    {
        let key = reverse_key(key);
        self.k_smallest_by_key(k, key)
    }
}

impl<T: ?Sized> ItertoolsExt for T where T: Iterator {}

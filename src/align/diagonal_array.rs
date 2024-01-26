/// A square diagonal array that is implemented as a single continuous vector
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiagonalArray<T> {
    len: usize,
    data: Vec<T>,
}

impl<T> DiagonalArray<T> {
    const fn length(n: usize) -> usize {
        (n + 1) * n / 2
    }

    pub unsafe fn get_unchecked(&self, index: [usize; 2]) -> &T {
        let index = Self::length(index[0]) + index[1];
        self.data.get_unchecked(index)
    }

    pub unsafe fn get_unchecked_mut(&mut self, index: [usize; 2]) -> &mut T {
        let index = Self::length(index[0]) + index[1];
        self.data.get_unchecked_mut(index)
    }
}

impl<T: Default + Clone> DiagonalArray<T> {
    /// Create a new diagonal array of the correct size, with all values initialised to the default value of the type
    pub fn new(len: usize) -> Self {
        Self {
            len,
            data: vec![T::default(); Self::length(len)],
        }
    }
}

impl<T> std::ops::Index<[usize; 2]> for DiagonalArray<T> {
    type Output = T;
    /// Index into the diagonal array
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        assert!(
            index[0] < self.len,
            "First index {} is outside of diagonal array with length {}",
            index[0],
            self.len
        );
        assert!(
            index[1] <= index[0],
            "Second index {} is outside of diagonal array with length {} at first index {}",
            index[1],
            self.len,
            index[0],
        );
        let index = Self::length(index[0]) + index[1];
        &self.data[index]
    }
}

impl<T> std::ops::IndexMut<[usize; 2]> for DiagonalArray<T> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        assert!(
            index[0] < self.len,
            "First index {} is outside of diagonal array with length {}",
            index[0],
            self.len
        );
        assert!(
            index[1] <= index[0],
            "Second index {} is outside of diagonal array with length {} at first index {}",
            index[1],
            self.len,
            index[0],
        );
        let index = Self::length(index[0]) + index[1];
        &mut self.data[index]
    }
}

#[cfg(test)]
mod tests {
    use super::DiagonalArray;

    #[test]
    fn create() {
        let mut array = DiagonalArray::<i8>::new(2);
        array[[0, 0]] = 1;
        array[[1, 0]] = 2;
        array[[1, 1]] = 3;
        assert_eq!(array[[0, 0]], 1);
        assert_eq!(array[[1, 0]], 2);
        assert_eq!(array[[1, 1]], 3);
    }
}

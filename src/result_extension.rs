pub trait ResultExtensions<T, E> {
    fn flat_err(self) -> Result<T, E>;
}

impl<T, E> ResultExtensions<T, E> for Result<T, Result<T, E>> {
    fn flat_err(self) -> Result<T, E> {
        match self {
            Ok(o) => Ok(o),
            Err(r) => r,
        }
    }
}

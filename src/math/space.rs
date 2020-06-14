use std::process::Output;

/// A trait representing any type that has a norm
pub trait NormedSpace {
    /// The associated scalar type: for a given type, there should only be one
    /// `impl` block corresponding to this trait, which is why it is not generic
    /// over `T`
    type Scalar;

    fn norm() -> Self::Scalar;
    fn norm_squared() -> Self::Scalar;
    fn normalize() -> Self;
}

/// A trait representing any type that has a norm AND defines a contraction inner product
///
/// See: http://www.jaapsuter.com/geometric-algebra.pdf
pub trait InnerSpace<T = Self>: NormedSpace {
    type Output;

    /// The dot (or "inner") product between an object of type `Self` and `T` (which defaults to `Self`)
    fn dot(self, rhs: T) -> Output;
}

pub trait OuterSpace<T = Self> {
    type Output;

    /// The wedge (or "outer") product between an object of type `Self` and `T` (which defaults to `Self`)
    fn wedge(self, rhs: T) -> Output;
}

/// A trait representing any multivector, which, for any other generic type `T`, must implement:
/// - The contraction inner product with elements of type `T`
/// - The outer product with elements of type `T`
/// - The geometric product with elements of type `T`: note that both the inner and outer products
///   can be defined in terms of the geometric product
pub trait Multivector<T = Self>: InnerSpace<T> + OuterSpace<T> {
    type Output;

    fn geometric_product(self, rhs: T) -> Output;
}

//https://github.com/rustgd/cgmath/blob/master/src/structure.rs
//https://crates.parity.io/alga/linear/trait.InnerSpace.html

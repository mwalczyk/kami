/// A trait representing any type that defines a contraction inner product (with
/// itself or another type).
pub trait Contraction<T = Self> {
    type Output;

    /// The dot (or "inner") product between an object of type `Self` and `T` (which defaults to `Self`).
    fn dot(self, rhs: T) -> Self::Output;
}

/// A trait representing any type that defines a exterior ("wedge") product (with
/// itself or another type).
pub trait Exterior<T = Self> {
    type Output;

    /// The wedge (or "exterior" or "outer") product between an object of type `Self` and `T`
    /// (which defaults to `Self`)
    fn wedge(self, rhs: T) -> Self::Output;
}

/// A trait representing any `k`-vector that has a norm. The norm is defined in such
/// a way that its square equals its contraction by its reversal. For example, the
/// squared norm of a bivector `a^b` would be given by `(b^a)•(a^b)`.
///
/// Given this definition, any `k`-vector that has a norm must also define a
/// contraction with objects of the same grade.
///
/// We separate these two traits, since a `k`-vector (i.e. blade) should only define
/// one norm, but it may define multiple versions of the contraction operation
/// (i.e. contraction between a bivector and a scalar, a bivector and a vector,
/// a bivector and another bivector, etc.).
///
/// Reference: "Understanding Geometric Algebra" (section 5.4)
pub trait Norm: Contraction<Self>
where
    Self: Sized,
{
    /// The associated scalar type: for a given type, there should only be one
    /// `impl` block corresponding to this trait, which is why it is not generic
    /// over `T`
    type Scalar;

    /// A measurement of the "magnitude" of this blade.
    fn norm(self) -> Self::Scalar;

    /// A measurement of the squared "magnitude" of this blade.
    fn norm_squared(self) -> Self::Scalar;

    /// Returns a "normalized" version of this blade, i.e. a new blade whose norm is 1.
    fn normalize(self) -> Self;
}

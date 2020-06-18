/// A trait representing any type that defines a contraction inner product (with
/// itself or another type).
///
/// Contracting a subspace by a `k`-vector lowers the dimensionality by `k`. Note
/// that subspaces of dimension less than 0 do not exist, so `(u^v)•α` (where `α`
/// is a scalar) equals 0. It follows that contraction of a scalar by *anything*
/// other than another scalar is 0 (non-existent). Similarly, contraction of a
/// vector by a bivector is 0, etc. Contraction of a `k`-vector with another
/// `k`-vector always results in a scalar that gives us some notion of similarity
/// (orientation and magnitude) between the two `k`-vectors.
///
/// Contraction of a subspace by another subspace is 0 if and only if the two
/// subspaces are orthogonal.
///
/// Reference: "Understanding Geometric Algebra" (section 5.3)
pub trait Contraction<T = Self> {
    type Output;

    /// The dot (or "inner") product between an object of type `Self` and `T`.
    fn dot(self, rhs: T) -> Self::Output;
}

/// A trait representing any type that defines a exterior ("wedge") product (with
/// itself or another type).
pub trait Exterior<T = Self> {
    type Output;

    /// The wedge (or "exterior" or "outer") product between an object of type `Self` and `T`.
    fn wedge(self, rhs: T) -> Self::Output;
}

/// A trait representing any `k`-vector that has a norm.
///
/// The norm of a `k`-vector can be defined as the contraction of a `k`-vector with
/// its reversal (i.e. `(v^u)•(u^v)`, in the case of bivectors). Note that scalars
/// and vectors don't have a reversal. As such, the norm of a scalar is simply the
/// square of that scalar. The norm of a vector is simply the dot product of that
/// vector with itself.
///
/// Expanding these definitions algebraically, we find that the norm of a bivector
/// is equal to the area of a parallelogram, the norm of a trivector is equal to
/// the volume of a parallelepiped, etc.
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

    /// Returns the inverse of this blade under the geometric product, such that `self`
    /// times its inverse is equal to 1.
    ///
    /// Reference: "Understanding Geometric Algebra" (section 6.3)
    fn inverse(self) -> Self;
}

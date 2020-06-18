use crate::math::space::{Exterior, Norm};

/// An enum representing the orientation of a `k`-blade. This has different
/// meanings depending on the grade of the blade. For example:
///
///     0-blades (scalars): or "oriented magnitude": positive or negative (sign)
///     1-blades (vectors) or "oriented length": positive or negative (direction)
///     2-blades (bivectors) or "oriented areas": clockwise or counterclockwise
///     3-blades (trivectors) or "oriented volumes": right-handed or left-handed
pub enum Orientation {
    Positive,
    Negative,
}

/// A trait representing any `k`-blade, where `k` refers to the dimension of the
/// subspace that the blade spans. For example:
///
///     0-blades -> scalars
///     1-blades -> vectors
///     2-blades -> trivectors
///     etc.
///
/// A `k`-blade is any object that can be expressed as the outer (exterior) product
/// of `k` vectors and is of grade `k`.
pub trait Blade: Norm + Exterior<Self> {
    /// The relative position / orientation of this blade, relative to some frame
    /// of reference (in this case, the Cartesian coordinate system). Currently
    /// unused.
    fn attitude() {
        unimplemented!()
    }

    /// Returns the orientation of this blade (positive or negative) which is the sign
    /// of its norm.
    fn orientation(&self) -> Orientation;

    /// Static function that returns the grade of this blade
    fn grade() -> usize;
}

/// Returns the total number of basis blades in a `dim`-dimensional space.
/// For example, a basis of Cl3 is:
///
///     { 1,    e1, e2, e3,    e12, e13, e23,    e123 }
///
/// which includes a basis scalar, 3 basis vectors, 3 basis bivectors, and
/// a basis trivector. The proof relies on combinatorial mathematics.
pub fn total_number_of_basis_blades(dim: usize) -> usize {
    2usize.pow(dim as u32)
}

/// Compute `n` factorial.
fn factorial(n: usize) -> usize {
    let mut fact = n;
    for i in 2..n {
        fact *= i;
    }
    return fact;
}

/// Compute `n` choose `k`, i.e. the binomial coefficient corresponding to
/// the parameters `n` and `k`.
fn binomial_coefficient(n: usize, k: usize) -> usize {
    let mut num = n;
    for i in 1..k {
        num *= n - i;
    }
    return num / factorial(k);
}

/// Returns the number of basis `k`-blades needed in an `n`-dimensional
/// space to represent arbitrary `k`-blades. For example, the number
/// of basis bivectors (2-blades) in 3-dimensional space is `number_of_basis_blades(3, 2)`,
/// which returns 3.
pub fn number_of_basis_blades(n: usize, k: usize) -> usize {
    binomial_coefficient(n, k)
}

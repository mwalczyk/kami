use std::ops::Add;

use num_traits::{One, Zero};

/// An object representing an oriented volume in 3-space. Trivectors form 3-dimensional
/// subspaces and are also known as the pseudo-scalars (of G3).
///
/// Trivectors are 3-blades. Each trivector can be decomposed onto a set of basis 3-blades.
/// In a 3-dimensional space, there is only 1 such basis 3-blade. However, in 4-dimensions,
/// for example, there are 4.
pub struct Trivector3<T> {
    /// 1st (and only) component, corresponding to the basis trivector `x^y^z`, i.e. e₁₂₃
    t1: T,
}

impl<T> Trivector3<T> {
    pub fn new(t1: T) -> Self {
        Self { t1 }
    }
}

impl<T> Trivector3<T>
where
    T: One,
{
    pub fn unit_xyz() -> Self {
        Trivector3::new(One::one())
    }
}

/// Add two trivectors, resulting in another trivector.
impl<T> Add for Trivector3<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Trivector3 {
            t1: self.t1 + other.t1,
        }
    }
}

/// Zero trivector.
impl<T: Zero> Zero for Trivector3<T> {
    #[inline]
    fn zero() -> Self {
        Trivector3 { t1: Zero::zero() }
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.t1.is_zero()
    }
}

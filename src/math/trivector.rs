use std::ops::Add;

use num_traits::{One, Zero};

/// An object representing an oriented volume in 3-space. Trivectors form 3-dimensional
/// subspaces and are also known as the pseudo-scalars of G3, since the vector space
/// formed by the trivectors is 1-dimensional. This follows from the fact that in G3, there
/// is only 1 basis trivector.
///
/// Trivectors are 3-blades. Each trivector can be decomposed onto a set of basis 3-blades.
/// In a 3-dimensional space, there is only 1 such basis 3-blade. However, in 4-dimensions,
/// for example, there are 4.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
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

impl<T: approx::AbsDiffEq> approx::AbsDiffEq for Trivector3<T>
where
    T::Epsilon: Copy,
{
    type Epsilon = T::Epsilon;

    fn default_epsilon() -> T::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: T::Epsilon) -> bool {
        T::abs_diff_eq(&self.t1, &other.t1, epsilon)
    }
}

impl<T: approx::RelativeEq> approx::RelativeEq for Trivector3<T>
where
    T::Epsilon: Copy,
{
    fn default_max_relative() -> T::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: T::Epsilon, max_relative: T::Epsilon) -> bool {
        T::relative_eq(&self.t1, &other.t1, epsilon, max_relative)
    }
}

impl<T: approx::UlpsEq> approx::UlpsEq for Trivector3<T>
where
    T::Epsilon: Copy,
{
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: T::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.t1, &other.t1, epsilon, max_ulps)
    }
}

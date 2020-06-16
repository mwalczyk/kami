use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::math::blade::{Blade, Orientation};
use crate::math::rotor::Rotor3;
use crate::math::trivector::Trivector3;
use crate::math::vector::Vector3;

use num_traits::{Float, One, Zero};

/// An object representing an oriented area in 3-space. Bivectors form 2-dimensional
/// subspaces.
///
/// Bivectors are also known as 2-blades.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Bivector3<T> {
    /// 1st component, corresponding to basis bivector `x^y` i.e. e₁₂
    pub b1: T,

    /// 2nd component, corresponding to basis bivector `x^z` i.e. e₁₃
    pub b2: T,

    /// 1st component, corresponding to basis bivector `y^z` i.e. e₂₃
    pub b3: T,
}

impl<T> Bivector3<T> {
    pub fn new(b1: T, b2: T, b3: T) -> Self {
        Self { b1, b2, b3 }
    }

    /// Constructs a trivector from a bivector `u^v` and a third vector `w` in R3.
    /// Here, `self` is the bivector `u^v` and `rhs` is `w`.
    ///
    /// Note that the resulting trivector is expressed in terms of the (single) basis
    /// trivector `x^y^z`. This is often written as `e₁₂₃`. A different basis could
    /// be used, such as `{ e₁₂₃ }`. We would simply need to switch some of the signs
    /// below to reflect this change.
    ///
    /// Another way to calculate this product would be via the determinant of the
    /// following 3x3 matrix:
    ///
    ///                 `| u₁ u₂ u₃ |`
    ///                 `| v₁ v₂ v₃ |`
    ///                 `| w₁ w₂ w₃ |`
    #[inline]
    pub fn wedge(self, rhs: Vector3<T>) -> Trivector3<T>
    where
        T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
    {
        Trivector3::new(self.b3 * rhs.a1 - self.b2 * rhs.a2 + self.b1 * rhs.a3)
    }
}

impl<T> Bivector3<T>
where
    T: One + Zero,
{
    /// Constructs the unit bivector in the `xy`-plane.
    pub fn unit_xy() -> Self {
        Self::new(One::one(), Zero::zero(), Zero::zero())
    }

    /// Constructs the unit bivector in the `xz`-plane.
    pub fn unit_xz() -> Self {
        Self::new(Zero::zero(), One::one(), Zero::zero())
    }

    /// Constructs the unit bivector in the `yz`-plane.
    pub fn unit_yz() -> Self {
        Self::new(Zero::zero(), Zero::zero(), One::one())
    }
}

impl<T> Bivector3<T>
where
    T: Add<Output = T> + Copy + Mul<Output = T>,
{
    #[inline]
    pub fn norm_squared(self) -> T {
        self.b1 * self.b1 + self.b2 * self.b2 + self.b3 * self.b3
    }

    #[inline]
    pub fn norm(self) -> T
    where
        T: Float,
    {
        self.norm_squared().sqrt()
    }

    /// Returns a normalized version of the bivector.
    #[inline]
    pub fn normalize(self) -> Self
    where
        T: Float,
    {
        self / self.norm()
    }
}

/// Add two bivectors, resulting in another bivector.
impl<T> Add for Bivector3<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    /// Two bivectors `B₁` and `B₂` can simply be added component-wise
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.b1 + rhs.b1, self.b2 + rhs.b2, self.b3 + rhs.b3)
    }
}

/// Divide a bivector by a scalar, resulting in another scaled bivector.
impl<T> Div<T> for Bivector3<T>
where
    T: Div<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self::new(self.b1 / rhs, self.b2 / rhs, self.b3 / rhs)
    }
}

/// Multiplying two bivectors, resulting in a rotor (i.e. a multivector with scalar
/// and bivector parts).
// TODO: is this correct? Does it make sense?
impl<T> Mul for Bivector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Neg<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Rotor3<T>;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        // The scalar part
        let scalar = -(self.b1 * rhs.b1 + self.b2 * rhs.b2 + self.b3 * rhs.b3);

        // The bivector part
        let bivector = Bivector3::new(
            self.b3 * rhs.b2 - self.b2 * rhs.b3,
            self.b1 * rhs.b3 - self.b3 * rhs.b1,
            self.b2 * rhs.b1 - self.b1 * rhs.b2,
        );

        Rotor3::new(scalar, bivector)
    }
}

/// Multiply a bivector by a scalar, resulting in another scaled bivector.
impl<T> Mul<T> for Bivector3<T>
where
    T: Mul<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Self::new(self.b1 * rhs, self.b2 * rhs, self.b3 * rhs)
    }
}

/// Negate a bivector, resulting in another bivector that is oriented in the opposite
/// direction.
impl<T> Neg for Bivector3<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.b1, -self.b2, -self.b3)
    }
}

/// Subtract two bivectors, resulting in another bivector.
impl<T> Sub for Bivector3<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    /// Two bivectors `B₁` and `B₂` can simply be subtracted component-wise
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.b1 - rhs.b1, self.b2 - rhs.b2, self.b3 - rhs.b3)
    }
}

/// Zero bivector.
impl<T> Zero for Bivector3<T>
where
    T: Zero,
{
    #[inline]
    fn zero() -> Self {
        Self::new(Zero::zero(), Zero::zero(), Zero::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.b1.is_zero() && self.b2.is_zero() && self.b3.is_zero()
    }
}

impl<T: approx::AbsDiffEq> approx::AbsDiffEq for Bivector3<T>
where
    T::Epsilon: Copy,
{
    type Epsilon = T::Epsilon;

    fn default_epsilon() -> T::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: T::Epsilon) -> bool {
        T::abs_diff_eq(&self.b1, &other.b1, epsilon)
            && T::abs_diff_eq(&self.b2, &other.b2, epsilon)
            && T::abs_diff_eq(&self.b3, &other.b3, epsilon)
    }
}

impl<T: approx::RelativeEq> approx::RelativeEq for Bivector3<T>
where
    T::Epsilon: Copy,
{
    fn default_max_relative() -> T::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: T::Epsilon, max_relative: T::Epsilon) -> bool {
        T::relative_eq(&self.b1, &other.b1, epsilon, max_relative)
            && T::relative_eq(&self.b2, &other.b2, epsilon, max_relative)
            && T::relative_eq(&self.b3, &other.b3, epsilon, max_relative)
    }
}

impl<T: approx::UlpsEq> approx::UlpsEq for Bivector3<T>
where
    T::Epsilon: Copy,
{
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: T::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.b1, &other.b1, epsilon, max_ulps)
            && T::ulps_eq(&self.b2, &other.b2, epsilon, max_ulps)
            && T::ulps_eq(&self.b3, &other.b3, epsilon, max_ulps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bivector_vector_product() {
        //let lhs = a * B;
        //let rhs = B * a;
        //assert_eq!(lhs, -rhs);
    }
}

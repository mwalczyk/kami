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
    /// 1st component, corresponding to basis bivector `x ^ y` i.e. e₁₂
    pub b1: T,

    /// 2nd component, corresponding to basis bivector `x ^ z` i.e. e₁₃
    pub b2: T,

    /// 1st component, corresponding to basis bivector `y ^ z` i.e. e₂₃
    pub b3: T,
}

impl<T> Bivector3<T> {
    pub fn new(b1: T, b2: T, b3: T) -> Self {
        Self { b1, b2, b3 }
    }

    /// Constructs a trivector from a bivector `u ^ v` and a third vector `w` in R3.
    /// Here, `self` is the bivector `u ^ v` and `rhs` is `w`.
    ///
    /// Note that the resulting trivector is expressed in terms of the (single) basis
    /// trivector `x ^ y ^ z`. This is often written as e_123. A different basis could
    /// be used, such as { e_132 }. We would simply need to switch some of the signs
    /// below to reflect this change.
    ///
    /// Another way to calculate this product would be via the determinant of the
    /// following 3x3 matrix:
    ///
    ///                 | u₁ u₂ u₃ |
    ///                 | v₁ v₂ v₃ |
    ///                 | w₁ w₂ w₃ |
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
    pub fn unit_xy() -> Self {
        Bivector3::new(One::one(), Zero::zero(), Zero::zero())
    }

    pub fn unit_xz() -> Self {
        Bivector3::new(Zero::zero(), One::one(), Zero::zero())
    }

    pub fn unit_yz() -> Self {
        Bivector3::new(Zero::zero(), Zero::zero(), One::one())
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

/// Add two bivectors, resulting in another bivector
impl<T: Add<Output = T>> Add for Bivector3<T> {
    type Output = Self;

    /// Two bivectors `B₁` and `B₂` can simply be added component-wise
    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Bivector3 {
            b1: self.b1 + other.b1,
            b2: self.b2 + other.b2,
            b3: self.b3 + other.b3,
        }
    }
}

/// Divide a bivector by a scalar, resulting in another scaled bivector
impl<T> Div<T> for Bivector3<T>
where
    T: Div<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Bivector3 {
            b1: self.b1 / rhs,
            b2: self.b2 / rhs,
            b3: self.b3 / rhs,
        }
    }
}

/// Multiply a bivector by a scalar, resulting in another scaled bivector
impl<T> Mul<T> for Bivector3<T>
where
    T: Mul<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Bivector3 {
            b1: self.b1 * rhs,
            b2: self.b2 * rhs,
            b3: self.b3 * rhs,
        }
    }
}

/// Negate a bivector, resulting in another bivector that is oriented in the opposite
/// direction
impl<T> Neg for Bivector3<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            b1: -self.b1,
            b2: -self.b2,
            b3: -self.b3,
        }
    }
}

/// Subtract two bivectors, resulting in another bivector
impl<T> Sub for Bivector3<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    /// Two bivectors `B₁` and `B₂` can simply be subtracted component-wise
    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Bivector3 {
            b1: self.b1 - other.b1,
            b2: self.b2 - other.b2,
            b3: self.b3 - other.b3,
        }
    }
}

/// Zero bivector
impl<T: Zero> Zero for Bivector3<T> {
    #[inline]
    fn zero() -> Self {
        Bivector3::new(Zero::zero(), Zero::zero(), Zero::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.b1.is_zero() && self.b2.is_zero() && self.b3.is_zero()
    }
}

/// Multiplying two bivectors, resulting in a rotor (i.e. a multivector with scalar
/// and bivector parts)
// TODO: is this correct? Does it make sense?
impl<T> Mul for Bivector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Neg<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Rotor3<T>;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        // The scalar part: notice the similarity to the dot product
        let a = -(self.b1 * other.b1 + self.b2 * other.b2 + self.b3 * other.b3);

        let b = Bivector3::new(
            self.b3 * other.b2 - self.b2 * other.b3,
            self.b1 * other.b3 - self.b3 * other.b1,
            self.b2 * other.b1 - self.b1 * other.b2,
        );

        Rotor3::new(a, b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero() {
        let a: Bivector3<f64> = Bivector3::zero();
        let b = Bivector3::new(0.0, 0.0, 0.0);
        assert_eq!(a, b);

        let a: Bivector3<i32> = Bivector3::zero();
        let b = Bivector3::new(0, 0, 0);
        assert_eq!(a, b);
    }

    fn bivector_vector_product() {
        //let lhs = a * B;
        //let rhs = B * a;
        //assert_eq!(lhs, -rhs);
    }
}

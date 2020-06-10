use std::ops::{Add, Div, Mul, Sub};

use crate::math::blade::{Blade, Orientation};
use crate::math::trivector::Trivector3;
use crate::math::vector::Vector3;

use num_traits::{One, Zero};

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
    fn wedge(self, rhs: Vector3<T>) -> Trivector3<T>
    where
        T: Mul<Output = T> + Add<Output = T> + Sub<Output = T> + Copy,
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

// Bivector-bivector multiplication: two bivectors can be multiplied, producing
// a multivector with scalar and bivector parts.
// impl<T: Mul<Output = T> + Copy> Mul for Bivector3<T> {
//     type Output = (T, Self);
//
//     #[inline]
//     fn mul(self, other: Self) -> Self::Output {
//         let scalar_part = Zero::zero();
//
//         let bivector_part = Bivector3 {
//             b_xy: self.b_xy * other.b_xy,
//             b_xz: self.b_xz * other.b_xz,
//             b_yz: self.b_yz * other.b_yz,
//         };
//
//         (scalar_part, bivector_part)
//     }
// }
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

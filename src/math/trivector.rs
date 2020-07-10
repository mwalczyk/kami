use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::math::space::{Contraction, Exterior, Norm};

use num_traits::{Float, One, Zero};

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
    /// 1st (and only) component, corresponding to the basis trivector , `e₁₂₃`, i.e. `x^y^z`
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
    /// Constructs the unit trivector.
    pub fn unit_xyz() -> Self {
        Trivector3::new(One::one())
    }
}

impl<T> Contraction for Trivector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Neg<Output = T>,
{
    type Output = T;

    /// Geometrically, we can interpret the dot product between two trivectors `(x^y^z)•(a^b^c)`
    /// as computing the intersection of the space `(a^b^c)` with the origin, which is orthogonal
    /// to the space `(x^y^z)`.
    ///
    /// Reference: "Understanding Geometric Algebra" (section 5.3)
    #[inline]
    fn dot(self, rhs: Self) -> Self::Output {
        -self.t1 * rhs.t1
    }
}

impl<T> Exterior for Trivector3<T>
where
    T: Zero,
{
    type Output = T;

    /// This is mostly here for completeness, as this will probably never be used. However, a
    /// trivector wedged with another trivector is always 0.
    #[inline]
    fn wedge(self, rhs: Self) -> Self::Output {
        T::zero()
    }
}

impl<T> Norm for Trivector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Neg<Output = T> + Copy + Float,
{
    type Scalar = T;

    /// Compute the (absolute value of the) volume of the trivector.
    #[inline]
    fn norm(self) -> Self::Scalar {
        self.norm_squared().sqrt()
    }

    /// For a trivector `T` formed by the wedge product `u^v^w`, the squared volume of `T`
    /// is given by `(w^v^u)•(u^v^w)`, i.e. contraction of the trivector by its reversal.
    #[inline]
    fn norm_squared(self) -> Self::Scalar {
        (-self).dot(self)
    }

    /// Returns a normalized version of the trivector, i.e. a trivector with the
    /// same orientation and attitude but unit volume.
    #[inline]
    fn normalize(self) -> Self {
        self / self.norm()
    }

    /// Returns the inverse `T⁻¹` of the trivector `T` under the geometric product
    /// such that `T⁻¹T` equals the scalar value 1.
    #[inline]
    fn inverse(self) -> Self {
        -self / self.norm_squared()
    }
}

impl<T> Add for Trivector3<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    /// Add two trivectors, resulting in another trivector.
    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Trivector3 {
            t1: self.t1 + other.t1,
        }
    }
}

impl<T> AddAssign for Trivector3<T>
where
    T: AddAssign,
{
    /// Add a trivector to this trivector.
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.t1 += rhs.t1;
    }
}

impl<T> Div<T> for Trivector3<T>
where
    T: Div<Output = T> + Copy,
{
    type Output = Self;

    /// Divide a trivector by a scalar, resulting in another scaled trivector.
    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self::new(self.t1 / rhs)
    }
}

impl<T> DivAssign<T> for Trivector3<T>
where
    T: DivAssign + Copy,
{
    /// Divide this trivector by a scalar.
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        self.t1 /= rhs;
    }
}

impl<T> Mul for Trivector3<T>
where
    T: Mul<Output = T> + Neg<Output = T>,
{
    type Output = T;

    /// Multiplying two trivectors, resulting in a scalar. The negative sign comes
    /// from a simple expansion of the geometric product.
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        -self.t1 * rhs.t1
    }
}

impl<T> MulAssign<T> for Trivector3<T>
where
    T: MulAssign + Copy,
{
    /// Multiply this trivector by a scalar.
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        self.t1 *= rhs;
    }
}

impl<T> Neg for Trivector3<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    /// Negate a trivector, resulting in another trivector that is oriented in the opposite
    /// direction.
    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.t1)
    }
}

impl<T> Sub for Trivector3<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    /// Subtract two trivectors, resulting in another trivector. Two trivectors `T₁`
    /// and `T₂` can simply be subtracted component-wise
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.t1 - rhs.t1)
    }
}

impl<T> SubAssign for Trivector3<T>
where
    T: SubAssign,
{
    /// Subtract a trivector from this trivector.
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.t1 -= rhs.t1;
    }
}

impl<T: Zero> Zero for Trivector3<T> {
    /// Zero trivector.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::vector::Vector3;

    #[test]
    fn norm() {
        let u = Vector3::new(1.0, 0.0, 0.0);
        let v = Vector3::new(0.0, 2.0, 0.0);
        let w = Vector3::new(0.0, 0.0, 3.0);
        let t = u.wedge(v).wedge(w);
        let res = t.norm();
        println!("Testing trivector norm...");
        println!("\t||T|| = {:?}", res);
        approx::assert_ulps_eq!(res, 6.0);
    }

    #[test]
    fn normalize() {
        let t = Trivector3::new(10.0);
        let res = t.normalize();
        println!("Testing trivector normalize with a positively oriented trivector...");
        println!("\tNormalized: {:?}", res);
        approx::assert_ulps_eq!(res, Trivector3::unit_xyz());

        let t = Trivector3::new(-10.0);
        let res = t.normalize();
        println!("Testing trivector normalize with a negatively oriented trivector...");
        println!("\tNormalized: {:?}", res);
        approx::assert_ulps_eq!(res, -Trivector3::unit_xyz());
    }

    #[test]
    fn inverse() {
        let t = Trivector3::new(10.0);
        let t_inv = t.inverse();
        let res_0 = t_inv * t;
        let res_1 = t * t_inv;
        println!("Testing trivector inverse...");
        println!("\tT⁻¹T = {:?}", res_0);
        println!("\tTT⁻¹ = {:?}", res_1);
        approx::assert_ulps_eq!(res_0, 1.0);
        approx::assert_ulps_eq!(res_1, 1.0);
    }

    #[test]
    fn geometric_product() {
        let t1 = Trivector3::new(10.0);
        let t2 = Trivector3::new(20.0);
        let res = t1 * t2;
        println!("Testing trivector geometric product...");
        println!("\tT₁T₂ = {:?}", res);
        approx::assert_ulps_eq!(res, -200.0);
    }
}

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::math::rotor::Rotor3;
use crate::math::space::{Contraction, Exterior, Norm};
use crate::math::trivector::Trivector3;
use crate::math::vector::Vector3;

use num_traits::{Float, One, Zero};

/// An object representing an oriented area in 3-space. Bivectors form 2-dimensional
/// subspaces.
///
/// Bivectors are also known as 2-blades.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Bivector3<T> {
    /// 1st component, corresponding to basis bivector `e₁₂`, i.e. `x^y`
    pub b1: T,

    /// 2nd component, corresponding to basis bivector `e₁₃`, i.e. `x^z`
    pub b2: T,

    /// 3rd component, corresponding to basis bivector `e₂₃`, i.e. `y^z`
    pub b3: T,
}

impl<T> Bivector3<T> {
    /// Constructs a new bivector from components.
    pub fn new(b1: T, b2: T, b3: T) -> Self {
        Self { b1, b2, b3 }
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

impl<T> Contraction for Bivector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Neg<Output = T>,
{
    type Output = T;

    /// Geometrically, we can interpret the dot product between two bivectors `(x^y)•(a^b)`
    /// as computing the intersection of the plane `(a^b)` with the surface normal to the
    /// plane `x^y`. The intersection is "strong" if they meet nearly perpendicularly and
    /// "weak" if they are nearly parallel.
    ///
    /// Reference: "Understanding Geometric Algebra" (section 5.3)
    #[inline]
    fn dot(self, rhs: Self) -> Self::Output {
        -(self.b1 * rhs.b1 + self.b2 * rhs.b2 + self.b3 * rhs.b3)
    }
}

impl<T> Exterior for Bivector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn wedge(self, rhs: Self) -> Self::Output {
        Bivector3::new(
            self.b3 * rhs.b2 - self.b2 * rhs.b3,
            self.b1 * rhs.b3 - self.b3 * rhs.b1,
            self.b2 * rhs.b1 - self.b1 * rhs.b2,
        )
    }
}

impl<T> Exterior<Vector3<T>> for Bivector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Trivector3<T>;

    /// Constructs a trivector from a bivector `u^v` and a third vector `w` in R3.
    /// Here, `self` is the bivector `u^v` and `rhs` is `w`.
    ///
    /// Note that the resulting trivector is expressed in terms of the (single) basis
    /// trivector `x^y^z`. This is often written as `e₁₂₃`. A different basis could
    /// be used, such as `{ e₁₃₂ }`. We would simply need to switch some of the signs
    /// below to reflect this change.
    ///
    /// Another way to calculate this product would be via the determinant of the
    /// following 3x3 matrix:
    ///
    ///                 `| u₁ u₂ u₃ |`
    ///                 `| v₁ v₂ v₃ |`
    ///                 `| w₁ w₂ w₃ |`
    #[inline]
    fn wedge(self, rhs: Vector3<T>) -> Self::Output {
        Trivector3::new(self.b3 * rhs.a1 - self.b2 * rhs.a2 + self.b1 * rhs.a3)
    }
}

impl<T> Norm for Bivector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Neg<Output = T> + Copy + Float,
{
    type Scalar = T;

    /// Compute the (absolute value of the) area of the bivector.
    #[inline]
    fn norm(self) -> Self::Scalar {
        self.norm_squared().sqrt()
    }

    /// For a bivector `B` formed by the wedge product `u^v`, the squared area of `B`
    /// is given by `(v^u)•(u^v)`, i.e. contraction of the bivector by its reversal.
    #[inline]
    fn norm_squared(self) -> Self::Scalar {
        (-self).dot(self)
    }

    /// Returns a normalized version of the bivector, i.e. a bivector with the
    /// same orientation and attitude but unit area.
    #[inline]
    fn normalize(self) -> Self {
        self / self.norm()
    }

    /// Returns the inverse `B⁻¹` of the bivector `B` under the geometric product
    /// such that `B⁻¹B` equals the scalar value 1.
    fn inverse(self) -> Self {
        -self / self.norm_squared()
    }
}

impl<T> Add for Bivector3<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    /// Add two bivectors, resulting in another bivector. Two bivectors `B₁`
    /// and `B₂` can simply be added component-wise
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.b1 + rhs.b1, self.b2 + rhs.b2, self.b3 + rhs.b3)
    }
}

impl<T> AddAssign for Bivector3<T>
where
    T: AddAssign,
{
    /// Add a bivector to this bivector.
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.b1 += rhs.b1;
        self.b2 += rhs.b2;
        self.b3 += rhs.b3;
    }
}

impl<T> Div<T> for Bivector3<T>
where
    T: Div<Output = T> + Copy,
{
    type Output = Self;

    /// Divide a bivector by a scalar, resulting in another scaled bivector.
    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self::new(self.b1 / rhs, self.b2 / rhs, self.b3 / rhs)
    }
}

impl<T> DivAssign<T> for Bivector3<T>
where
    T: DivAssign + Copy,
{
    /// Divide this bivector by a scalar.
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        self.b1 /= rhs;
        self.b2 /= rhs;
        self.b3 /= rhs;
    }
}

impl<T> Mul for Bivector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Neg<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Rotor3<T>;

    /// Multiplying two bivectors, resulting in a rotor (i.e. a multivector with scalar
    /// and bivector parts).
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

impl<T> Mul<T> for Bivector3<T>
where
    T: Mul<Output = T> + Copy,
{
    type Output = Self;

    /// Multiply a bivector by a scalar, resulting in another scaled bivector.
    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Self::new(self.b1 * rhs, self.b2 * rhs, self.b3 * rhs)
    }
}

impl<T> MulAssign<T> for Bivector3<T>
where
    T: MulAssign + Copy,
{
    /// Multiply this bivector by a scalar.
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        self.b1 *= rhs;
        self.b2 *= rhs;
        self.b3 *= rhs;
    }
}

impl<T> Neg for Bivector3<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    /// Negate a bivector, resulting in another bivector that is oriented in the opposite
    /// direction. The negation (or reversal) of a bivector formed by `u^v` is simply
    /// `v^u`.
    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.b1, -self.b2, -self.b3)
    }
}

impl<T> Sub for Bivector3<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    /// Subtract two bivectors, resulting in another bivector. Two bivectors `B₁`
    /// and `B₂` can simply be subtracted component-wise
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.b1 - rhs.b1, self.b2 - rhs.b2, self.b3 - rhs.b3)
    }
}

impl<T> SubAssign for Bivector3<T>
where
    T: SubAssign,
{
    /// Subtract a bivector from this bivector.
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.b1 -= rhs.b1;
        self.b2 -= rhs.b2;
        self.b3 -= rhs.b3;
    }
}

impl<T> Zero for Bivector3<T>
where
    T: Zero,
{
    /// Zero bivector.
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

    #[test]
    fn dot() {
        let xy = Bivector3::<f32>::unit_xy();
        let xz = Bivector3::unit_xz();
        let res = xy.dot(xz);
        println!("Testing bivector dot product with orthogonal planes...");
        println!("\t||(x^y)•(x^z)|| = {:?}", res);
        approx::assert_ulps_eq!(res, 0.0);

        // Note that when two bivectors are parallel, the dot product returns -1, not 1
        // (as is the case with two parallel vectors)
        let res = xy.dot(xy);
        println!("Testing bivector dot product with parallel planes...");
        println!("\t||(x^y)•(x^y)|| = {:?}", res);
        approx::assert_ulps_eq!(res, -1.0);
    }

    #[test]
    fn norm() {
        // First, some basic sanity checks
        let xy = Bivector3::<f32>::unit_xy();
        let xz = Bivector3::<f32>::unit_xz();
        let yz = Bivector3::<f32>::unit_yz();
        approx::assert_ulps_eq!(xy.norm(), 1.0);
        approx::assert_ulps_eq!(xz.norm(), 1.0);
        approx::assert_ulps_eq!(yz.norm(), 1.0);

        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res = u.wedge(v);

        // These should all be equal:
        // 1) The norm of `B`
        // 2) The (algebraic expansion of the) contraction of `B` with its reversal
        // 3) The contraction (dot product) of `B` with its reversal
        // 4) The scalar part of the geometric product between `B` and its reversal
        //
        // Reference: Reference: "Understanding Geometric Algebra" (section 5.4)
        let norm = res.norm_squared();
        let via_contraction = u.norm_squared() * v.norm_squared() - u.dot(v).powf(2.0);
        let via_dot_product = (v.wedge(u)).dot(u.wedge(v));
        let via_geometric_product = v.wedge(u) * u.wedge(v);
        println!("Testing bivector norm...");
        println!("\tCalculated (via norm function): {:?}, ", norm);
        println!("\tVia contraction: {:?}", via_contraction);
        println!("\tVia dot product: {:?}", via_dot_product);
        println!(
            "\tVia geometric product (scalar part) <(v^u)(u^v)>₀ = {:?}",
            via_geometric_product.scalar
        );
        approx::assert_ulps_eq!(norm, via_contraction);
        approx::assert_ulps_eq!(norm, via_dot_product);
        approx::assert_ulps_eq!(norm, via_geometric_product.scalar);
    }

    #[test]
    fn normalize() {
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res_0 = u.wedge(v).normalize().norm();
        let res_1 = v.wedge(u).normalize().norm();
        println!("Testing bivector normalize...");
        println!("\t||u^v|| = {:?}", res_0);
        println!("\t||v^u|| = {:?}", res_1);
        approx::assert_ulps_eq!(res_0, 1.0);
        approx::assert_ulps_eq!(res_1, 1.0);
    }

    #[test]
    fn inverse() {
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let b = u.wedge(v);
        let b_inv = b.inverse();
        let res_0 = b_inv * b;
        let res_1 = b * b_inv;
        println!("Testing bivector inverse...");
        println!("\tB⁻¹B = {:?}", res_0);
        println!("\tBB⁻¹ = {:?}", res_1);
        approx::assert_ulps_eq!(res_0.scalar, 1.0);
        approx::assert_ulps_eq!(res_0.bivector, Bivector3::zero());
        approx::assert_ulps_eq!(res_1.scalar, 1.0);
        approx::assert_ulps_eq!(res_1.bivector, Bivector3::zero());
    }

    #[test]
    fn geometric_product() {
        // Construct the first bivector
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res_0 = u.wedge(v);

        // Construct the second bivector
        let a = Vector3::new(6.0, 7.0, 8.0);
        let b = Vector3::new(9.0, 10.0, 11.0);
        let res_1 = a.wedge(b);

        // Compute the geometric product between the two bivectors directly
        let res = res_0 * res_1;

        // Compute the geometric product using the alternate definition as the sum
        // of the dot product and wedge product between the two bivectors
        let from_parts = Rotor3::new(res_0.dot(res_1), res_0.wedge(res_1));

        println!("Testing bivector geometric product with another bivector...");
        println!("\t(u^v)(a^b) = {:?}", res);
        println!("\tFrom dot and wedge products: {:?}", from_parts);
        approx::assert_ulps_eq!(res, from_parts);
    }
}

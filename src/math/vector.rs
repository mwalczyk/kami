use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::math::bivector::Bivector3;
use crate::math::rotor::Rotor3;
use crate::math::space::{Contraction, Exterior, Norm};
use crate::math::trivector::Trivector3;

use num_traits::{Float, One, Zero};

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Vector3<T> {
    /// 1st component, corresponding to basis vector `e₁`, i.e. the positive `x`-axis
    pub a1: T,

    /// 2nd component, corresponding to basis vector `e₂`, i.e. the positive `y`-axis
    pub a2: T,

    /// 3rd component, corresponding to basis vector `e₃`, i.e. the positive `z`-axis
    pub a3: T,
}

impl<T> Vector3<T> {
    /// Construct a new vector with the specified components.
    pub fn new(a1: T, a2: T, a3: T) -> Self {
        Self { a1, a2, a3 }
    }
}

impl<T> Vector3<T>
where
    T: One + Zero,
{
    /// Construct a unit vector aligned with the positive `x`-axis.
    pub fn unit_x() -> Self {
        Vector3::new(One::one(), Zero::zero(), Zero::zero())
    }

    /// Construct a unit vector aligned with the positive `y`-axis.
    pub fn unit_y() -> Self {
        Vector3::new(Zero::zero(), One::one(), Zero::zero())
    }

    /// Construct a unit vector aligned with the positive `z`-axis.
    pub fn unit_z() -> Self {
        Vector3::new(Zero::zero(), Zero::zero(), One::one())
    }
}

impl<T> Contraction for Vector3<T>
where
    T: Add<Output = T> + Mul<Output = T>,
{
    type Output = T;

    /// Also known as the "inner product," the dot product is a measure of
    /// similarity between two vectors. It takes as input two vectors `u`
    /// and `v` and returns a scalar, i.e. a grade-0 object. This is the
    /// contraction of a vector with itself.
    #[inline]
    fn dot(self, rhs: Self) -> Self::Output {
        self.a1 * rhs.a1 + self.a2 * rhs.a2 + self.a3 * rhs.a3
    }
}

impl<T> Exterior for Vector3<T>
where
    T: Mul<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Bivector3<T>;

    /// Constructs a bivector from a pair of vectors in R3 via the outer
    /// (wedge) product: `u^v`.
    ///
    /// Note that the resulting bivector is expressed in terms of the
    /// three basis bivectors `x^y`, `x^z`, and `y^z`. These are
    /// often written as `e₁₂`, `e₁₃`, and `e₂₃`, respectively. A different
    /// basis could be used, such as `{ e₁₂, e₂₃, e₃₁ }`. We would
    /// simply need to switch some of the signs below to reflect this
    /// change.
    ///
    /// Another way to calculate this product is via the (pseudo)
    /// determinant of the following 3x3 matrix:
    ///
    ///             `| e₂₃  e₁₃  e₁₂ |`
    ///             `|  u₁   u₂   u₃ |`
    ///             `|  v₁   v₂   v₃ |`
    #[inline]
    fn wedge(self, rhs: Self) -> Self::Output {
        Bivector3::new(
            self.a1 * rhs.a2 - self.a2 * rhs.a1, // XY (i.e. e₁₂)
            self.a1 * rhs.a3 - self.a3 * rhs.a1, // XZ (i.e. e₁₃)
            self.a2 * rhs.a3 - self.a3 * rhs.a2, // YZ (i.e. e₂₃)
        )
    }
}

impl<T> Norm for Vector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Copy + Float,
{
    type Scalar = T;

    /// Compute the (absolute value of the) length of the vector `u`.
    #[inline]
    fn norm(self) -> Self::Scalar {
        self.norm_squared().sqrt()
    }

    /// For a vector `u`, the squared length of `u` is given by `u•u`.
    #[inline]
    fn norm_squared(self) -> Self::Scalar {
        self.dot(self)
    }

    /// Returns a normalized version of the vector, i.e. a vector with the
    /// same orientation and attitude but unit length.
    #[inline]
    fn normalize(self) -> Self {
        self / self.norm()
    }

    /// The inverse of the vector under the geometric product.
    #[inline]
    fn inverse(self) -> Self {
        self / self.norm_squared()
    }
}

impl<T> Vector3<T>
where
    T: Float,
{
    /// Calculates the angle between two vectors. If `u` and `v` are unit
    /// vectors, then the angle between them is simply `cosθ = u•v`.
    #[inline]
    pub fn angle(self, rhs: Self) -> T {
        let u_norm = self.normalize();
        let v_norm = rhs.normalize();
        u_norm.dot(v_norm).acos()
    }
}

impl<T> Vector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Debug + Copy + Float,
{
    /// Returns the projection of `u` onto `v`, using the dot product
    /// and the geometric inverse of `v`.
    #[inline]
    pub fn project(self, rhs: Self) -> Self {
        // Normally, we would write this as `(u•v)v⁻¹`, but we can't multiply
        // a scalar by a vector on the left with the current API
        rhs.inverse() * self.dot(rhs)
    }

    /// Returns the reflection of `u` about `v`. The geometric algebra version of
    /// this problem involves 2 geometric products: `u' = v⁻¹uv`. In general,
    /// this sequence of operations will either produce a vector or a trivector.
    /// However, the trivector case only occurs when the 3 vectors are linearly
    /// independent, which is obviously not the case with the reflection formula
    /// listed above. See the footnotes of https://marctenbosch.com/quaternions/
    /// for details.
    pub fn reflect(self, rhs: Self) -> Self
    where
        T: Div<Output = T> + Add<Output = T> + Copy + Mul<Output = T> + Sub<Output = T>,
    {
        // The trivector part should always be 0
        let rotor = rhs.inverse() * self;
        let (vector, trivector) = rotor * rhs;
        //println!("{:?}", trivector);

        vector
    }
}

impl<T> Add for Vector3<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    /// Add two vectors, resulting in another vector.
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.a1 + rhs.a1, self.a2 + rhs.a2, self.a3 + rhs.a3)
    }
}

impl<T> AddAssign for Vector3<T>
where
    T: AddAssign,
{
    /// Add a vector to this vector.
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.a1 += rhs.a1;
        self.a2 += rhs.a2;
        self.a3 += rhs.a3;
    }
}

impl<T> Div<T> for Vector3<T>
where
    T: Div<Output = T> + Copy,
{
    type Output = Self;

    /// Divide a vector by a scalar, resulting in another scaled vector.
    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self::new(self.a1 / rhs, self.a2 / rhs, self.a3 / rhs)
    }
}

impl<T> DivAssign<T> for Vector3<T>
where
    T: DivAssign + Copy,
{
    /// Divide this vector by a scalar.
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        self.a1 /= rhs;
        self.a2 /= rhs;
        self.a3 /= rhs;
    }
}

impl<T> Mul for Vector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Rotor3<T>;

    /// Multiply two vectors via the geometric product, producing a multivector
    /// with a scalar (grade-0) part and a bivector (grade-2) part, a.k.a a rotor.
    fn mul(self, rhs: Self) -> Self::Output {
        Rotor3::new(self.dot(rhs), self.wedge(rhs))
    }
}

impl<T> Mul<T> for Vector3<T>
where
    T: Mul<Output = T> + Copy,
{
    type Output = Self;

    /// Multiply a vector by a scalar, resulting in another scaled vector.
    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Self::new(self.a1 * rhs, self.a2 * rhs, self.a3 * rhs)
    }
}

impl<T> MulAssign<T> for Vector3<T>
where
    T: MulAssign + Copy,
{
    /// Multiply this vector by a scalar.
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        self.a1 *= rhs;
        self.a2 *= rhs;
        self.a3 *= rhs;
    }
}

impl<T> Mul<Bivector3<T>> for Vector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Neg<Output = T> + Sub<Output = T> + Copy,
{
    type Output = (Vector3<T>, Trivector3<T>);

    /// Multiply a vector and a bivector via the geometric product, resulting in
    /// a multivector with a grade-1 (vector) part `<aB>₁` and a grade-3 (trivector)
    /// part `<aB>₃`.
    #[inline]
    fn mul(self, rhs: Bivector3<T>) -> Self::Output {
        // The vector (grade_1) part
        let vector = Vector3::new(
            -self.a2 * rhs.b1 - self.a3 * rhs.b2,
            self.a1 * rhs.b1 - self.a3 * rhs.b3,
            self.a1 * rhs.b2 + self.a2 * rhs.b3,
        );

        // The trivector (grade-3) part
        let trivector = Trivector3::new(self.a1 * rhs.b3 - self.a2 * rhs.b2 + self.a3 * rhs.b1);

        (vector, trivector)
    }
}

impl<T> Neg for Vector3<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    /// Negate a vector, resulting in another vector that is oriented in the opposite
    /// direction.
    fn neg(self) -> Self::Output {
        Self::new(-self.a1, -self.a2, -self.a3)
    }
}

impl<T> Sub for Vector3<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    /// Subtract two vectors, resulting in another vector.
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.a1 - rhs.a1, self.a2 - rhs.a2, self.a3 - rhs.a3)
    }
}

impl<T> SubAssign for Vector3<T>
where
    T: SubAssign,
{
    /// Subtract a vector from this vector.
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.a1 -= rhs.a1;
        self.a2 -= rhs.a2;
        self.a3 -= rhs.a3;
    }
}

impl<T> Zero for Vector3<T>
where
    T: Zero,
{
    /// Zero vector.
    #[inline]
    fn zero() -> Self {
        Self::new(Zero::zero(), Zero::zero(), Zero::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.a1.is_zero() && self.a2.is_zero() && self.a3.is_zero()
    }
}

impl<T: approx::AbsDiffEq> approx::AbsDiffEq for Vector3<T>
where
    T::Epsilon: Copy,
{
    type Epsilon = T::Epsilon;

    fn default_epsilon() -> T::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: T::Epsilon) -> bool {
        T::abs_diff_eq(&self.a1, &other.a1, epsilon)
            && T::abs_diff_eq(&self.a2, &other.a2, epsilon)
            && T::abs_diff_eq(&self.a3, &other.a3, epsilon)
    }
}

impl<T: approx::RelativeEq> approx::RelativeEq for Vector3<T>
where
    T::Epsilon: Copy,
{
    fn default_max_relative() -> T::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: T::Epsilon, max_relative: T::Epsilon) -> bool {
        T::relative_eq(&self.a1, &other.a1, epsilon, max_relative)
            && T::relative_eq(&self.a2, &other.a2, epsilon, max_relative)
            && T::relative_eq(&self.a3, &other.a3, epsilon, max_relative)
    }
}

impl<T: approx::UlpsEq> approx::UlpsEq for Vector3<T>
where
    T::Epsilon: Copy,
{
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: T::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.a1, &other.a1, epsilon, max_ulps)
            && T::ulps_eq(&self.a2, &other.a2, epsilon, max_ulps)
            && T::ulps_eq(&self.a3, &other.a3, epsilon, max_ulps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot() {
        // Take the dot product between the two vectors `u` and `v`
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res = u.dot(v);

        // The inner (dot) product and outer (wedge) product are subsidiary
        // operations, related to the geometric product as follows:
        //
        // u•v = ½(uv + vu)
        // u^v = ½(uv - vu)
        //
        // These are rotors that we can add component-wise
        let uv = u * v;
        let vu = v * u;
        let rotor = (uv + vu) * 0.5;
        println!("Testing vector dot product and its relation to the geometric product...");
        println!("\tu•v = {:?}", res);
        println!("\t½(uv + vu) = {:?}", rotor);
        approx::assert_ulps_eq!(res, 14.0);
        approx::assert_ulps_eq!(res, rotor.scalar);
        approx::assert_ulps_eq!(rotor.bivector, Bivector3::zero());
    }

    #[test]
    fn wedge() {
        // `u` and `v` are the same, so `u^v` should return the zero bivector
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(0.0, 1.0, 2.0);
        let res_0 = u.wedge(v);
        let res_1 = v.wedge(u);
        println!("Testing vector wedge product with parallel vectors...");
        println!("\tu^v = {:?}", res_0);
        println!("\tv^u = {:?}", res_1);
        approx::assert_ulps_eq!(res_0, Bivector3::zero());
        approx::assert_ulps_eq!(res_1, Bivector3::zero());

        // `u` and `v` are orthonormal to each other, so `u^v` should return a unit
        // bivector in one of the basis planes
        let u = Vector3::new(1.0, 0.0, 0.0);
        let v = Vector3::new(0.0, 1.0, 0.0);
        let res_0 = u.wedge(v);
        let res_1 = v.wedge(u);
        println!("Testing vector wedge product with orthogonal vectors...");
        println!("\tu^v = {:?}", res_0);
        println!("\tv^u = {:?}", res_1);
        approx::assert_ulps_eq!(res_0, Bivector3::new(1.0, 0.0, 0.0));
        approx::assert_ulps_eq!(res_1, Bivector3::new(-1.0, 0.0, 0.0));

        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res = u.wedge(v);

        // The inner (dot) product and outer (wedge) product are subsidiary
        // operations, related to the geometric product as follows:
        //
        // u•v = ½(uv + vu)
        // u^v = ½(uv - vu)
        //
        // These are rotors that we can add component-wise
        let uv = u * v;
        let vu = v * u;
        let rotor = (uv - vu) * 0.5;
        println!("Testing vector wedge product and its relation to the geometric product...");
        println!("\tu^v = {:?}", res);
        println!("\t½(uv - vu) = {:?}", rotor);
        approx::assert_ulps_eq!(0.0, rotor.scalar);
        approx::assert_ulps_eq!(res, rotor.bivector);
    }

    #[test]
    fn inverse() {
        // For any vector `u`, taking the geometric product of `u⁻¹` and `u`
        // should be 1 (a scalar with no bivector part)
        let u = Vector3::new(1.0, 2.0, 3.0);
        let u_inv = u.inverse();
        let res = u_inv * u;
        println!("Testing vector inverse...");
        println!("\tu⁻¹u = {:?}", res);
        approx::assert_ulps_eq!(res.scalar, 1.0);
        approx::assert_ulps_eq!(res.bivector, Bivector3::zero());
    }

    #[test]
    fn angle() {
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res = u.angle(v);
        println!("Testing vector angle...");
        println!("\tAngle between u and v: {:?}", res);
        approx::assert_ulps_eq!(res, ((7.0 * 10.0.sqrt()) / 25.0).acos(), max_ulps = 8);
    }

    #[test]
    fn project() {
        let u = Vector3::new(1.0, 1.0, 0.0);
        let v = Vector3::new(2.0, 0.0, 0.0);
        let res = u.project(v);
        println!("Testing vector projection...");
        println!("\t{:?}", res);
        approx::assert_ulps_eq!(res, Vector3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn reflect() {
        let u = Vector3::<f32>::unit_x();
        let v = Vector3::unit_y();
        let res = u.reflect(v);
        println!("Testing vector reflection with basis vectors...");
        println!("\tReflected vector: {:?}", res);

        // Reflect `u` through `v` using the geometric algebra approach
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res = u.reflect(v);

        // Reflect `u` through `v` using the "traditional" linear algebra approach: `2(v•u)v - u`,
        // where the vector `v` *must* be normalized
        //
        // Reference: https://www.fabrizioduroni.it/2017/08/25/how-to-calculate-reflection-vector.html
        let v_normalized = v.normalize();
        let traditional = v_normalized * (2.0 * u.dot(v_normalized)) - u;
        println!("Testing vector reflection with arbitrary vectors...");
        println!("\tReflected vector: {:?}", res);
        println!("\tVia traditional formula (non-GA): {:?}", traditional);
        approx::assert_ulps_eq!(res, traditional, max_ulps = 6);
        approx::assert_ulps_eq!(res, Vector3::new(42.0 / 25.0, 31.0 / 25.0, 4.0 / 5.0));
    }

    #[test]
    fn add() {
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res = u + v;
        println!("Testing vector addition...");
        println!("\tu + v = {:?}", res);
        approx::assert_ulps_eq!(res, Vector3::new(3.0, 5.0, 7.0));
    }

    #[test]
    fn sub() {
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res = u - v;
        println!("Testing vector subtraction...");
        println!("\tu - v = {:?}", res);
        approx::assert_ulps_eq!(res, Vector3::new(-3.0, -3.0, -3.0));
    }
}

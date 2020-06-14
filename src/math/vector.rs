use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::math::bivector::Bivector3;
use crate::math::rotor::Rotor3;
use crate::math::trivector::Trivector3;

use num_traits::{Float, One, Zero};

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Vector3<T> {
    /// 1st component, corresponding to basis vector `x` i.e. e₁
    pub a1: T,

    /// 2nd component, corresponding to basis vector `y` i.e. e₂
    pub a2: T,

    /// 3rd component, corresponding to basis vector `z` i.e. e₃
    pub a3: T,
}

impl<T> Vector3<T> {
    pub fn new(a1: T, a2: T, a3: T) -> Self {
        Self { a1, a2, a3 }
    }
}

impl<T> Vector3<T>
where
    T: One + Zero,
{
    pub fn unit_x() -> Self {
        Vector3::new(One::one(), Zero::zero(), Zero::zero())
    }

    pub fn unit_y() -> Self {
        Vector3::new(Zero::zero(), One::one(), Zero::zero())
    }

    pub fn unit_z() -> Self {
        Vector3::new(Zero::zero(), Zero::zero(), One::one())
    }
}

impl<T> Vector3<T>
where
    T: Mul<Output = T> + Sub<Output = T> + Copy,
{
    /// Constructs a bivector from a pair of vectors in R3 via the outer
    /// (wedge) product: `u^v`. Here, `self` is `u` and `rhs` is `v`.
    ///
    /// Note that the resulting bivector is expressed in terms of the
    /// three basis bivectors `x^y`, `x^z`, and `y^z`. These are
    /// often written as e₁₂, e₁₃, and e₂₃, respectively. A different
    /// basis could be used, such as { e₁₂, e₂₃, e₃₁ }. We would
    /// simply need to switch some of the signs below to reflect this
    /// change.
    ///
    /// Another way to calculate this product would be via the (pseudo)
    /// determinant of the following 3x3 matrix:
    ///
    ///             `| e₂₃  e₁₃  e₁₂ |`
    ///             `| u₁   u₂   u₃  |`
    ///             `| v₁   v₂   v₃  |`
    #[inline]
    pub fn wedge(self, rhs: Vector3<T>) -> Bivector3<T> {
        Bivector3::new(
            self.a1 * rhs.a2 - self.a2 * rhs.a1, // XY (i.e. e₁₂)
            self.a1 * rhs.a3 - self.a3 * rhs.a1, // XZ (i.e. e₁₃)
            self.a2 * rhs.a3 - self.a3 * rhs.a2, // YZ (i.e. e₂₃)
        )
    }
}

impl<T> Vector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Copy,
{
    /// Also known as the "dot product," the inner product is a measure of
    /// similarity between two vectors. It takes as input two vectors `u`
    /// and `v` and returns a scalar, i.e. grade-0 object.
    #[inline]
    pub fn dot(self, rhs: Self) -> T {
        self.a1 * rhs.a1 + self.a2 * rhs.a2 + self.a3 * rhs.a3
    }

    /// For a vector `u`, the squared length of `u` is given by `u•u`.
    #[inline]
    pub fn norm_squared(self) -> T {
        self.dot(self)
    }

    /// Compute the length of a vector `u`: note that a square root
    /// operation is required, and thus, this function only works
    /// for floats.
    #[inline]
    pub fn norm(self) -> T
    where
        T: Float,
    {
        self.norm_squared().sqrt()
    }

    /// The inverse of the vector under the geometric product.
    #[inline]
    pub fn inverse(self) -> Self
    where
        T: Div<Output = T>,
    {
        self / self.norm_squared()
    }

    /// Returns a normalized version of the vector.
    #[inline]
    pub fn normalize(self) -> Self
    where
        T: Float,
    {
        self / self.norm()
    }

    /// Returns the projection of `u` onto `v`, using the dot product
    /// and the geometric inverse of `v`.
    #[inline]
    pub fn project(self, rhs: Self) -> Self
    where
        T: Div<Output = T>,
    {
        // Normally, we would write this as `(u•v)v⁻¹`, but we can't multiply
        // a scalar by a vector on the left with the current API
        rhs.inverse() * self.dot(rhs)
    }

    // Returns the reflection of `u` about `v`. The geometric algebra version of
    // this problem involves 2 geometric products: `u' = v⁻¹uv`. In general,
    // this sequence of operations will either produce a vector or a trivector.
    // However, the trivector case only occurs when the 3 vectors are linearly
    // independent, which is obviously not the case with the reflection formula
    // listed above. See the footnotes of https://marctenbosch.com/quaternions/
    // for details.
    // fn reflect(self, rhs: Self) -> Self
    // where
    //     T: Div<Output = T> + Add<Output = T> + Copy + Mul<Output = T> + Sub<Output = T>,
    // {
    //     let (scalar, bivector) = rhs.inverse() * self;
    // }
}

/// Add two vectors, resulting in another vector
impl<T> Add for Vector3<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            a1: self.a1 + rhs.a1,
            a2: self.a2 + rhs.a2,
            a3: self.a3 + rhs.a3,
        }
    }
}

/// Divide a vector by a scalar, resulting in another scaled vector
impl<T> Div<T> for Vector3<T>
where
    T: Div<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self {
            a1: self.a1 / rhs,
            a2: self.a2 / rhs,
            a3: self.a3 / rhs,
        }
    }
}

/// Multiply a vector by a scalar, resulting in another scaled vector
impl<T> Mul<T> for Vector3<T>
where
    T: Mul<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Self {
            a1: self.a1 * rhs,
            a2: self.a2 * rhs,
            a3: self.a3 * rhs,
        }
    }
}

/// Multiply two vectors via the geometric product, producing a multivector
/// with a scalar (grade-0) part and a bivector (grade-2) part, a.k.a a rotor
impl<T> Mul for Vector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Rotor3<T>;

    fn mul(self, rhs: Self) -> Self::Output {
        Rotor3::new(self.dot(rhs), self.wedge(rhs))
    }
}

/// Multiply a vector and a bivector via the geometric product, producing a
/// multivector with a vector (grade-1) part and a trivector (grade-3) part
impl<T> Mul<Bivector3<T>> for Vector3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Neg<Output = T> + Sub<Output = T> + Copy,
{
    type Output = (Vector3<T>, Trivector3<T>);

    /// Multiplying a vector `a` and a bivector `B` results in a multivector with
    /// a grade-1 (vector) part `<aB>₁` and a grade-3 (trivector) part `<aB>₃`.
    #[inline]
    fn mul(self, rhs: Bivector3<T>) -> Self::Output {
        // The vector (grade_1) part
        let scalar = Vector3::new(
            -self.a2 * rhs.b1 - self.a3 * rhs.b2,
            self.a1 * rhs.b1 - self.a3 * rhs.b3,
            self.a1 * rhs.b2 + self.a2 * rhs.b3,
        );

        // The trivector (grade-3) part
        let trivector = Trivector3::new(self.a1 * rhs.b3 - self.a2 * rhs.b2 + self.a3 * rhs.b1);

        (scalar, trivector)
    }
}

/// Negate a vector, resulting in another vector that is oriented in the opposite
/// direction
impl<T> Neg for Vector3<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            a1: -self.a1,
            a2: -self.a2,
            a3: -self.a3,
        }
    }
}

/// Subtract two vectors, resulting in another vector
impl<T> Sub for Vector3<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            a1: self.a1 - rhs.a1,
            a2: self.a2 - rhs.a2,
            a3: self.a3 - rhs.a3,
        }
    }
}

/// Zero vector
impl<T> Zero for Vector3<T>
where
    T: Zero,
{
    #[inline]
    fn zero() -> Self {
        Self::new(Zero::zero(), Zero::zero(), Zero::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.a1.is_zero() && self.a2.is_zero() && self.a3.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add() {
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res = u + v;
        println!("Testing vector addition...");
        println!("  u + v = {:?}", res);
        assert_eq!(res, Vector3::new(3.0, 5.0, 7.0));
    }

    #[test]
    fn sub() {
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(3.0, 4.0, 5.0);
        let res = u - v;
        println!("Testing vector subtraction...");
        println!("  u - v = {:?}", res);
        assert_eq!(res, Vector3::new(-3.0, -3.0, -3.0));
    }

    #[test]
    fn outer_product() {
        // `u` and `v` are the same, so `u^v` should return the zero bivector
        let u = Vector3::new(0.0, 1.0, 2.0);
        let v = Vector3::new(0.0, 1.0, 2.0);
        let res_0 = u.wedge(v);
        let res_1 = v.wedge(u);
        println!("Testing vector outer product with parallel vectors...");
        println!("  u^v = {:?}", res_0);
        println!("  v^u = {:?}", res_1);
        assert_eq!(res_0, Bivector3::zero());
        assert_eq!(res_1, Bivector3::zero());

        let u = Vector3::new(1.0, 0.0, 0.0);
        let v = Vector3::new(0.0, 1.0, 0.0);
        let res_0 = u.wedge(v);
        let res_1 = v.wedge(u);
        println!("Testing vector outer product with orthogonal vectors...");
        println!("  u^v = {:?}", res_0);
        println!("  v^u = {:?}", res_1);
        assert_eq!(res_0, Bivector3::new(1.0, 0.0, 0.0));
        assert_eq!(res_1, Bivector3::new(-1.0, 0.0, 0.0));

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
        let rotor = Rotor3::new(uv.scalar - vu.scalar, uv.bivector - vu.bivector) * 0.5;
        println!("Testing vector outer product and its relation to the geometric product...");
        println!("  u^v = {:?}", res);
        println!("  ½(uv - vu) = {:?}", rotor);
        assert_eq!(0.0, rotor.scalar);
        assert_eq!(res, rotor.bivector);
    }

    #[test]
    fn inner_product() {
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
        let rotor = Rotor3::new(uv.scalar + vu.scalar, uv.bivector + vu.bivector) * 0.5;
        println!("Testing vector inner product and its relation to the geometric product...");
        println!("  u•v = {:?}", res);
        println!("  ½(uv + vu) = {:?}", rotor);
        assert!((res - 14.0).abs() <= 0.001);
        assert!((rotor.scalar - res).abs() <= 0.001);
        assert_eq!(rotor.bivector, Bivector3::zero());
    }

    #[test]
    fn inverse() {
        // For any vector `u`, taking the geometric product of `u⁻¹` and `u`
        // should be 1 (a scalar with no bivector part)
        let u = Vector3::new(1.0, 2.0, 3.0);
        let u_inv = u.inverse();
        let res = u_inv * u;
        println!("Testing vector inverse...");
        println!("  u⁻¹u = {:?}", res);
        assert!((res.scalar - 1.0).abs() <= 0.001);
        assert_eq!(res.bivector, Bivector3::zero());
    }

    #[test]
    fn project() {
        let u = Vector3::new(1.0, 1.0, 0.0);
        let v = Vector3::new(2.0, 0.0, 0.0);
        let res = u.project(v);
        println!("Testing vector projection...");
        println!("  {:?}", res);
        assert_eq!(res, Vector3::new(1.0, 0.0, 0.0));
    }
}

use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::math::bivector::Bivector3;
use crate::math::rotor::Rotor3;
use crate::math::trivector::Trivector3;

use num_traits::{Float, Zero};

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
    fn new(a1: T, a2: T, a3: T) -> Self {
        Self { a1, a2, a3 }
    }
}

impl<T> Vector3<T>
where
    T: Copy + Mul<Output = T> + Sub<Output = T>,
{
    /// Constructs a bivector from a pair of vectors in R3 via the outer
    /// (wedge) product: `u ^ v`. Here, `self` is `u` and `rhs` is `v`.
    ///
    /// Note that the resulting bivector is expressed in terms of the
    /// three basis bivectors `x ^ y`, `x ^ z`, and `y ^ z`. These are
    /// often written as e₁₂, e₁₃, and e₂₃, respectively. A different
    /// basis could be used, such as { e₁₂, e₂₃, e₃₁ }. We would
    /// simply need to switch some of the signs below to reflect this
    /// change.
    ///
    /// Another way to calculate this product would be via the (pseudo)
    /// determinant of the following 3x3 matrix:
    ///
    ///             | e₂₃  e₁₃  e₁₂ |
    ///             | u₁   u₂   u₃  |
    ///             | v₁   v₂   v₃  |
    #[inline]
    fn wedge(self, rhs: Vector3<T>) -> Bivector3<T> {
        Bivector3::new(
            self.a1 * rhs.a2 - self.a2 * rhs.a1, // XY (i.e. e₁₂)
            self.a1 * rhs.a3 - self.a3 * rhs.a1, // XZ (i.e. e₁₃)
            self.a2 * rhs.a3 - self.a3 * rhs.a2, // YZ (i.e. e₂₃)
        )
    }
}

impl<T> Vector3<T>
where
    T: Add<Output = T> + Copy + Mul<Output = T>,
{
    /// Also known as the "dot product," the inner product is a measure of
    /// similarity between two vectors. It takes as input two vectors `u`
    /// and `v` and returns a scalar, i.e. grade-0 object.
    #[inline]
    fn dot(self, rhs: Self) -> T {
        self.a1 * rhs.a1 + self.a2 * rhs.a2 + self.a3 * rhs.a3
    }

    /// For a vector `u`, the squared length of `u` is given by `u.u`.
    #[inline]
    fn norm_squared(self) -> T {
        self.dot(self)
    }

    /// Compute the length of a vector `u`: note that a square root
    /// operation is required, and thus, this function only works
    /// for floats.
    #[inline]
    fn norm(self) -> T
    where
        T: Float,
    {
        self.norm_squared().sqrt()
    }

    /// The inverse of the vector under the geometric product.
    #[inline]
    fn inverse(self) -> Self
    where
        T: Div<Output = T>,
    {
        self / self.norm_squared()
    }

    /// Returns a normalized version of the vector.
    #[inline]
    fn normalize(self) -> Self
    where
        T: Float,
    {
        self / self.norm()
    }

    /// Returns the projection of `u` onto `v`, using the dot product
    /// and the geometric inverse of `v`.
    #[inline]
    fn project(self, rhs: Self) -> Self
    where
        T: Div<Output = T>,
    {
        // Normally, we would write this as `(u.v)v_inv`, but we can't multiply
        // a scalar by a vector on the left
        rhs.inverse() * self.dot(rhs)
    }

    // Returns the reflection of `u` about `v`. The geometric algebra version of
    // this problem involves 2 geometric products: `u' = v_inv * u * v`. In general,
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
        Vector3 {
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
        Vector3 {
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
        Vector3 {
            a1: self.a1 * rhs,
            a2: self.a2 * rhs,
            a3: self.a3 * rhs,
        }
    }
}

/// Multiply two vectors via the geometric product, producing a multivector
/// with a scalar (grade-0) part and a bivector (grade-2) part, a.k.a a rotor
impl<T> Mul<Vector3<T>> for Vector3<T>
where
    T: Add<Output = T> + Copy + Mul<Output = T> + Sub<Output = T>,
{
    type Output = Rotor3<T>;

    fn mul(self, rhs: Vector3<T>) -> Self::Output {
        Rotor3::new(self.dot(rhs), self.wedge(rhs))
    }
}

/// Multiply a vector and a bivector via the geometric product, producing a
/// multivector with a vector (grade-1) part and a trivector (grade-3) part
impl<T> Mul<Bivector3<T>> for Vector3<T>
where
    T: Add<Output = T> + Copy + Mul<Output = T> + Neg<Output = T> + Sub<Output = T>,
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
        Vector3 {
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
        Vector3::new(Zero::zero(), Zero::zero(), Zero::zero())
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
        let a = Vector3::new(0.0, 1.0, 2.0);
        let b = Vector3::new(3.0, 4.0, 5.0);
        let res = a + b;
        println!("Testing vector addition...");
        assert_eq!(res, Vector3::new(3.0, 5.0, 7.0));
    }

    #[test]
    fn sub() {
        let a = Vector3::new(0.0, 1.0, 2.0);
        let b = Vector3::new(3.0, 4.0, 5.0);
        let res = a - b;
        println!("Testing vector subtraction...");
        assert_eq!(res, Vector3::new(-3.0, -3.0, -3.0));
    }

    #[test]
    fn outer_product() {
        // `a` and `b` are the same, so `a ^ b` should return the zero bivector
        let a = Vector3::new(0.0, 1.0, 2.0);
        let b = Vector3::new(0.0, 1.0, 2.0);
        let res_0 = a.wedge(b);
        let res_1 = b.wedge(a);
        println!("Testing vector outer product with parallel vectors...");
        println!("  a^b = {:?}", res_0);
        println!("  b^a = {:?}", res_1);
        assert_eq!(res_0, Bivector3::zero());
        assert_eq!(res_1, Bivector3::zero());

        let a = Vector3::new(1.0, 0.0, 0.0);
        let b = Vector3::new(0.0, 1.0, 0.0);
        let res_0 = a.wedge(b);
        let res_1 = b.wedge(a);
        println!("Testing vector outer product with basis vectors e1 and e2...");
        println!("  a^b = {:?}", res_0);
        println!("  b^a = {:?}", res_1);
        assert_eq!(res_0, Bivector3::new(1.0, 0.0, 0.0));
        assert_eq!(res_1, Bivector3::new(-1.0, 0.0, 0.0));
    }

    fn inner_product() {
        let a = Vector3::new(0.0, 1.0, 2.0);
        let b = Vector3::new(3.0, 4.0, 5.0);

        let res = a * b;

        // The inner (dot) product and outer (wedge) product are subsidiary
        // operations, related to the geometric product as follows:
        //
        // u.v = 1/2 (uv + vu)
        // u^v = 1/2 (uv - vu)
    }

    #[test]
    fn inverse() {
        // For any vector `u`, taking the geometric product of `u` inverse and `u`
        // should be 1 (a scalar with no bivector part)
        let a = Vector3::new(1.0f64, 2.0, 3.0);
        let b = a.inverse();

        let res = b * a;
        println!("Testing vector inverse...");
        println!("  b_inv*a = {:?}", res);
        assert!((res.a - 1.0).abs() <= 0.001);
        assert_eq!(res.b, Bivector3::zero());
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

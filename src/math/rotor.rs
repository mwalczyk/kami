use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::math::bivector::Bivector3;
use crate::math::space::{Contraction, Exterior};
use crate::math::vector::Vector3;

use num_traits::{Float, One, Zero};

/// A rotor is the composition of a scalar, grade-0 object and a bivector,
/// grade-2 object. Although scalars and bivectors can't be added together,
/// we can treat rotors the same way we treat complex numbers, where the
/// "real" and "imaginary" parts are kept together but cannot be further
/// simplified.
///
/// We call `uv = u.v + u^v` a rotor because if we multiply by `uv` on
/// both sides of a vector we perform a rotation.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Rotor3<T> {
    pub scalar: T,
    pub bivector: Bivector3<T>,
}

impl<T> Rotor3<T> {
    /// Construct a rotor from a scalar and a bivector.
    pub fn new(scalar: T, bivector: Bivector3<T>) -> Self {
        Self { scalar, bivector }
    }

    /// Construct the "identity" rotor, i.e. a rotor that does not induce
    /// any rotation.
    pub fn identity() -> Self
    where
        T: One + Zero,
    {
        Self::new(T::one(), Bivector3::zero())
    }

    /// Returns the conjugate of this rotor (i.e. its negation). This is a
    /// rotor that performs a rotation in the same plane but opposite direction.
    pub fn conjugate(self) -> Self
    where
        T: Neg<Output = T>,
    {
        Self::new(self.scalar, -self.bivector)
    }
}

impl<T> Rotor3<T>
where
    T: Float,
{
    /// Construct a rotor by taking the geometric product of two vectors `u`
    /// and `v`. The geometric product produces a rotor with a scalar (grade-0)
    /// part and a bivector (grade-2) part.
    ///
    /// The resulting rotor represents a rotation in the `uv`-plane by an
    /// angle that is equal to *twice* the angle of separation between
    /// `u` and `v`. In other words, it is the rotor that rotates from `u`
    /// to `v`.
    ///
    /// Usually, we want to construct a rotor that represents a rotation by
    /// angle `θ/2` - that is, the "true" angle between `u` and `v`. We can
    /// achieve this by adding 1 to the scalar part of the rotor, as shown
    /// below.
    pub fn from_vectors(u: Vector3<T>, v: Vector3<T>) -> Self
    where
        T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
    {
        // Note that we don't actually use the geometric product function here,
        // since we need to add 1 to the scalar part first
        //
        // Note also that we flip the order of `u` and `v` below: this is because
        // the formula for rotating a vector with a rotor is:
        //
        //          `x' = (vu)x(uv)`
        //
        // In most of the literature, the left-hand multiplier `vu` corresponds to
        // the rotor `R`, whereas the right-hand multiplier `uv` corresponds to the
        // inverse rotor, `R†` (we want to construct `R`)
        let scalar = T::one() + v.dot(u);
        let bivector = v.wedge(u);

        // Make sure to normalize
        Self::new(scalar, bivector).normalize()
    }

    /// Construct a rotor from an angle (in radians) and a plane of rotation (bivector).
    /// Note that the bivector must be normalized in order for this to work.
    ///
    /// This construction follows from Euler's formula:
    ///
    ///             `R = cosθ - sinθI`
    ///
    /// where `I` is the unit bivector.
    pub fn from_angle_plane(angle_radians: T, plane: Bivector3<T>) -> Self {
        let sin = (angle_radians / T::from(2.0).unwrap()).sin();
        let cos = (angle_radians / T::from(2.0).unwrap()).cos();
        Self::new(cos, plane * -sin)
    }

    /// Construct a rotor that performs a rotation of `angle_radians` in the `xy`-plane.
    pub fn from_xy_rotation(angle_radians: T) -> Self {
        Self::from_angle_plane(angle_radians, Bivector3::unit_xy())
    }

    /// Construct a rotor that performs a rotation of `angle_radians` in the `xz`-plane.
    pub fn from_xz_rotation(angle_radians: T) -> Self {
        Self::from_angle_plane(angle_radians, Bivector3::unit_xz())
    }

    /// Construct a rotor that performs a rotation of `angle_radians` in the `yz`-plane.
    pub fn from_yz_rotation(angle_radians: T) -> Self {
        Self::from_angle_plane(angle_radians, Bivector3::unit_yz())
    }
}

impl<T> Rotor3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Copy,
{
    /// The squared norm of this rotor.
    #[inline]
    pub fn norm_squared(self) -> T {
        self.scalar * self.scalar + self.bivector.norm_squared()
    }

    /// The norm of a rotor is simply the sum of the norms of its scalar and
    /// bivector parts.
    #[inline]
    pub fn norm(self) -> T
    where
        T: Float,
    {
        self.norm_squared().sqrt()
    }

    /// Returns a normalized version of the rotor.
    #[inline]
    fn normalize(self) -> Self
    where
        T: Float,
    {
        self / self.norm()
    }

    /// To rotate a vector using a rotor, we apply the formula `x' = RxR†`.
    pub fn rotate(self, u: Vector3<T>) -> Vector3<T>
    where
        T: Sub<Output = T> + Neg<Output = T>,
    {
        // First, do: x' = Rx
        let a = self.scalar;
        let Bivector3 { b1, b2, b3 } = self.bivector;
        let Vector3 {
            a1: u1,
            a2: u2,
            a3: u3,
        } = u;

        let q = Vector3::new(
            a * u1 + b1 * u2 + b2 * u3,
            a * u2 - b1 * u1 + b3 * u3,
            a * u3 - b2 * u1 - b3 * u2,
        );

        // Trivector part
        let q123 = -b1 * u3 + b2 * u2 - b3 * u1; // TODO: accounted for sign flip here

        // Then, do: x'' = x'R†
        let Vector3 {
            a1: q1,
            a2: q2,
            a3: q3,
        } = q;
        let y = q123; // TODO: had to add a negative sign here

        Vector3::new(
            a * q1 + q2 * b1 + q3 * b2 - y * b3, // TODO: had to swap the sign of the `y` term in this line
            a * q2 - q1 * b1 + y * b2 + q3 * b3, //   and this line
            a * q3 - y * b1 - q1 * b2 - q2 * b3,
        )
    }
}

/// Add two rotors, resulting in another rotor.
impl<T> Add for Rotor3<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    /// Two rotors `R₁` and `R₂` can simply be added component-wise
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.scalar + rhs.scalar, self.bivector + rhs.bivector)
    }
}

/// Multiply two rotors, resulting in another rotor that represents the combined rotation
/// of `R₁` and `R₂`. Note that `R₁R₂` implies that the rotation `R₂` happens first. The
/// calculations below are simply the expansion of the geometric product between the two
/// rotors.
impl<T> Mul for Rotor3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        // The scalar part
        let scalar = self.scalar * rhs.scalar
            - self.bivector.b1 * rhs.bivector.b1
            - self.bivector.b2 * rhs.bivector.b2
            - self.bivector.b3 * rhs.bivector.b3;

        // The bivector part
        let bivector = Bivector3::new(
            self.scalar * rhs.bivector.b1 + rhs.scalar * self.bivector.b1
                - self.bivector.b2 * rhs.bivector.b3
                + self.bivector.b3 * rhs.bivector.b2,
            self.scalar * rhs.bivector.b2
                + rhs.scalar * self.bivector.b2
                + self.bivector.b1 * rhs.bivector.b3
                - self.bivector.b3 * rhs.bivector.b1,
            self.scalar * rhs.bivector.b3 + rhs.scalar * self.bivector.b3
                - self.bivector.b1 * rhs.bivector.b2
                + self.bivector.b2 * rhs.bivector.b1,
        );

        Self::new(scalar, bivector)
    }
}

/// Multiply a rotor by a scalar, resulting in another scaled rotor.
impl<T> Mul<T> for Rotor3<T>
where
    T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Self::new(self.scalar * rhs, self.bivector * rhs)
    }
}

/// Divide a rotor by a scalar, resulting in another scaled rotor.
impl<T> Div<T> for Rotor3<T>
where
    T: Div<Output = T> + Copy,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self::new(self.scalar / rhs, self.bivector / rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity() {
        let x = Vector3::<f32>::unit_x();
        let rotor = Rotor3::identity();

        let x_prime = rotor.rotate(x);

        println!("Testing rotor identity...");
        println!("\tRotor: {:?}", rotor);
        println!("\tRotated vector: {:?}", x_prime);
        approx::assert_ulps_eq!(x_prime, Vector3::unit_x());
    }

    #[test]
    fn rotate() {
        // Input vectors
        let x_0 = Vector3::new(1.0, 0.0, 0.0);
        let x_1 = Vector3::new(1.0 / 2.0.sqrt(), 1.0 / 2.0.sqrt(), 0.0);
        let x_2 = Vector3::new(1.0 / 3.0.sqrt(), 1.0 / 3.0.sqrt(), 1.0 / 3.0.sqrt());

        // Construct a rotor in the xy-plane that performs a 90-degree rotation
        let u = Vector3::new(1.0, 0.0, 0.0);
        let v = Vector3::new(0.0, 1.0, 0.0);
        let rotor = Rotor3::from_vectors(u, v);

        let x_0_prime = rotor.rotate(x_0);
        let x_1_prime = rotor.rotate(x_1);
        let x_2_prime = rotor.rotate(x_2);

        println!("Testing rotor rotation (basic)...");
        println!("\tRotor: {:?}", rotor);
        println!("\tRotated vector #0: {:?}", x_0_prime);
        println!("\tRotated vector #1: {:?}", x_1_prime);
        println!("\tRotated vector #2: {:?}", x_2_prime);

        // Should be mapped to the positive y-axis
        approx::assert_ulps_eq!(x_0_prime, Vector3::new(0.0, 1.0, 0.0));

        // Should be mapped to the vector that runs diagonally (upper-left quadrant) of the xy-plane
        approx::assert_ulps_eq!(
            x_1_prime,
            Vector3::new(-1.0 / 2.0.sqrt(), 1.0 / 2.0.sqrt(), 0.0)
        );

        // Construct a rotor in the yz-plane that performs a 45-degree rotation
        let rotor = Rotor3::from_angle_plane(std::f32::consts::FRAC_PI_4, Bivector3::unit_yz());

        let x_0_prime = rotor.rotate(x_0);
        let x_1_prime = rotor.rotate(x_1);
        let x_2_prime = rotor.rotate(x_2);

        println!("Testing rotor rotation (angle + plane)...");
        println!("\tRotor: {:?}", rotor);
        println!("\tRotated vector #0: {:?}", x_0_prime);
        println!("\tRotated vector #1: {:?}", x_1_prime);
        println!("\tRotated vector #2: {:?}", x_2_prime);

        // Should be mapped to the positive x-axis (i.e. the rotation should have no effect)
        approx::assert_ulps_eq!(x_0_prime, Vector3::new(1.0, 0.0, 0.0));

        // Should be mapped to a vector whose x component is unchanged
        approx::assert_ulps_eq!(x_1_prime, Vector3::new(1.0 / 2.0.sqrt(), 0.5, 0.5));

        // Should be mapped to a vector whose x component is unchanged and whose y component is close to 0
        approx::assert_ulps_eq!(
            x_2_prime,
            Vector3::new(1.0 / 3.0.sqrt(), 0.0, (2.0 / 3.0).sqrt())
        );
    }
}

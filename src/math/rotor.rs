use crate::math::bivector::Bivector3;

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
    pub a: T,
    pub b: Bivector3<T>,
}

impl<T> Rotor3<T> {
    pub fn new(a: T, b: Bivector3<T>) -> Self {
        Self { a, b }
    }
}
// Construct a rotor from a scalar and a bivector

// Construct a rotor from two vectors `u` and `v` by taking the geometric
// product of `u` and `v`

// Multiply two rotors, resulting in another rotor

// Multiply a rotor by a vector, resulting in a rotated vector

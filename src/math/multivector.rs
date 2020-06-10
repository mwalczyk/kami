use crate::math::bivector::Bivector3;
use crate::math::trivector::Trivector3;
use crate::math::vector::Vector3;

struct Multivector<T> {
    scalar: T,
    vector: Vector3<T>,
    bivector: Bivector3<T>,
    trivector: Trivector3<T>,
}

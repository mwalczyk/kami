use crate::math::bivector::Bivector3;
use crate::math::space::Exterior;
use crate::math::vector::Vector3;

use approx::UlpsEq;
use num_traits::Zero;

pub struct Line {
    point: Vector3<f32>,
    direction: Vector3<f32>,
}

impl Line {
    /// Constructs a new (infinite) line that contains the point `point` and has slope `direction`.
    pub fn new(point: &Vector3<f32>, direction: &Vector3<f32>) -> Line {
        Line {
            point: *point,
            direction: *direction,
        }
    }

    /// Returns `true` if the specified point is on the line and `false` otherwise. A point
    /// `x` is on the line if `(x-a)^u = 0`. In other words, since `a` is on the line, `x`
    /// will be on the line if and only if the vector `x-a` (the vector from `a` to `x`) is
    /// parallel to the vector `u` (in which case, the bivector formed between them will have
    /// zero area).
    ///
    /// Reference: http://geocalc.clas.asu.edu/GA_Primer/GA_Primer/introduction-to-geometric/high-school-geometry-with/solutions-for-straight-line.html
    pub fn contains(&self, x: &Vector3<f32>) -> bool {
        let bivector = (*x - self.point).wedge(self.direction);
        bivector.ulps_eq(&Bivector3::zero(), std::f32::EPSILON, 4)
    }

    pub fn intersect(&self, rhs: &Line) -> Option<Vector3<f32>> {
        // a = (q^v) / (u^v)
        let a = 1.0;

        // b = (p^u) / (v^u);
        let b = 1.0;

        // x = au - bv

        None
    }
}

struct LineSegment {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contains() {
        // Create a line that passes through (1, 1, 0) with slope (1, 2, 0)
        let point = Vector3::new(1.0, 1.0, 0.0);
        let direction = Vector3::new(1.0, 2.0, 0.0);
        let line = Line::new(&point, &direction);
        let x_0 = point + direction;
        let x_1 = point - direction * 2.0;
        let x_2 = Vector3::new(-1.0, 1.0, 0.0);
        let res_0 = line.contains(&x_0);
        let res_1 = line.contains(&x_1);
        let res_2 = line.contains(&x_2);

        println!("Testing line contains in 2D...");
        println!("\tIs on line x₀: {:?}", res_0);
        println!("\tIs on line x₁: {:?}", res_1);
        println!("\tIs on line x₂: {:?}", res_2);
        assert_eq!(res_0, true);
        assert_eq!(res_1, true);
        assert_eq!(res_2, false);
    }
}

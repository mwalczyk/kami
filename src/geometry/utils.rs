use nalgebra_glm::Vec2;

const EPS: f32 = 0.001;

/// Determines whether a triangle is degenerate (i.e. has one or more sides with length 0)
pub fn is_degenerate_triangle(a: &Vec2, b: &Vec2, c: &Vec2) -> bool {
    let side_a = nalgebra_glm::distance(b, c);
    let side_b = nalgebra_glm::distance(a, c);
    let side_c = nalgebra_glm::distance(a, b);

    if side_a.abs() < EPS || side_b.abs() < EPS || side_c.abs() < EPS {
        return true;
    }

    false
}

/// Calculates the incenter of the triangle formed by three points `a`,
/// `b`, and `c`. Returns `None` if the triangle is degenerate.
pub fn calculate_triangle_incenter(a: &Vec2, b: &Vec2, c: &Vec2) -> Option<Vec2> {
    if is_degenerate_triangle(a, b, c) {
        return None;
    }

    let side_a = nalgebra_glm::distance(b, c);
    let side_b = nalgebra_glm::distance(a, c);
    let side_c = nalgebra_glm::distance(a, b);
    let inv_p = 1.0 / (side_a + side_b + side_c);

    let incenter = Vec2::new(
        (side_a * a.x + side_b * b.x + side_c * c.x) * inv_p,
        (side_a * a.y + side_b * b.y + side_c * c.y) * inv_p,
    );

    Some(incenter)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_degenerate_triangle() {}
}

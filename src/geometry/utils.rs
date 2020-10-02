use nalgebra_glm::Vec2;

/// Determines whether the vector from `a` to `b` is (nearly) parallel to the vector pointing
/// from `c` to `d` (within an epsilon).
pub fn parallel(a: &Vec2, b: &Vec2, c: &Vec2, d: &Vec2) -> bool {
    let dir_ab = b - a;
    let dir_cd = d - c;

    // If AB and CD are nearly parallel, this will be either -1 or 1, so we take
    // the absolute value
    let abs_dot = nalgebra_glm::normalize_dot(&dir_ab, &dir_cd).abs();

    if (abs_dot - 1.0).abs() < nalgebra_glm::epsilon() {
        return true;
    }

    false
}

/// Determines whether a triangle is degenerate (i.e. has one or more sides with length 0).
pub fn is_degenerate_triangle(a: &Vec2, b: &Vec2, c: &Vec2) -> bool {
    let side_a = nalgebra_glm::distance(b, c);
    let side_b = nalgebra_glm::distance(a, c);
    let side_c = nalgebra_glm::distance(a, b);

    // Triangle inequality: could also sort `{A, B, C}` and simply check if `(A + B) <= C`,
    // but this is more readable for now
    if (side_a + side_b) <= side_c || (side_a + side_c) <= side_b || (side_b + side_c) <= side_a {
        return true;
    }

    false
}

/// Calculates the incenter of the triangle formed by three points `a`,
/// `b`, and `c`. Returns `None` if the triangle is degenerate.
pub fn calculate_triangle_incenter(a: &Vec2, b: &Vec2, c: &Vec2) -> Option<(Vec2, f32)> {
    // First, check if this is a valid triangle
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

    // The semiperimeter of the triangle
    let p = (side_a + side_b + side_c) * 0.5;

    // Heron's formula for calculating the area of an arbitrary triangle
    let area = (p * (p - side_a) * (p - side_b) * (p - side_c)).sqrt();

    // Finally, calculate the inradius from the two quantities above
    let inradius = area / p;

    Some((incenter, inradius))
}

/// Calculates the intersection between two line segments (`ab` and `cd`) or `None` if they do not intersect.
pub fn calculate_line_segment_intersection(a: &Vec2, b: &Vec2, c: &Vec2, d: &Vec2) -> Option<Vec2> {
    let delta_x0 = b.x - a.x;
    let delta_y0 = b.y - a.y;
    let delta_x1 = d.x - c.x;
    let delta_y1 = d.y - c.y;
    let denom = (-delta_x1 * delta_y0 + delta_x0 * delta_y1);

    let s = (-delta_y0 * (a.x - c.x) + delta_x0 * (a.y - c.y)) / denom;
    let t = (delta_x1 * (a.y - c.y) - delta_y1 * (a.x - c.x)) / denom;

    if s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0 {
        return Some(Vec2::new(a.x + (t * delta_x0), a.y + (t * delta_y0)));
    }

    None
}

/// Given a line segment `ab` and a point `p` on the line, calculates a point `q` such that the segment
/// `pq` is perpendicular to `ab`.
pub fn calculate_perpendicular_point(a: &Vec2, b: &Vec2, p: &Vec2) -> Vec2 {
    let numer = (b.y - a.y) * (p.x - a.x) - (b.x - a.x) * (p.y - a.y);
    let denom = (b.y - a.y).powf(2.0) + (b.x - a.x).powf(2.0);
    let k = numer / denom;

    Vec2::new(p.x - k * (b.y - a.y), p.y + k * (b.x - a.x))
}

/// Calculates the midpoint along the segment `ab`.
pub fn calculate_midpoint(a: &Vec2, b: &Vec2) -> Vec2 {
    (*a + *b) * 0.5
}

/// Calculates the perpendicular bisector of segment `ab`. In particular, this function will return
/// a pair of points: the midpoint of the segment `ab`, and a point `offset`-units away from `ab` in
/// a direction perpendicular to `ab`.
pub fn calculate_perpendicular_bisector(a: &Vec2, b: &Vec2, offset: f32) -> (Vec2, Vec2) {
    // First, calculate the midpoint of the line segment from `a` to `b`
    let midpoint = calculate_midpoint(a, b);

    // Then, calculate the slope of the line segment and the slope of the perpendicular
    let slope = (b.y - a.y) / (b.x - a.x);
    let slope_perpendicular = -1.0 / slope;
    let intercept = midpoint.y - slope_perpendicular * midpoint.x;

    // Add a small offset in the x-direction (this is arbitrary)
    let offset_point = midpoint + Vec2::new(offset, 0.0);
    let second = Vec2::new(
        offset_point.x,
        offset_point.x * slope_perpendicular + intercept,
    );

    (midpoint, second)
}

/// Checks if the specified point `p` lies along the line segment `ab`.
pub fn is_on_line_segment(a: &Vec2, b: &Vec2, p: &Vec2, include_endpoints: bool) -> bool {
    // The total length of the line segment
    let ab = nalgebra_glm::distance(a, b);

    // The lengths of the two sub-segments joining each endpoint to the specified point
    let ap = nalgebra_glm::distance(a, p);
    let pb = nalgebra_glm::distance(p, b);

    // Return `false` if the specified point is one of the endpoints of the line segment
    ((ap + pb) - ab).abs() < nalgebra_glm::epsilon()
        && ap.abs() > nalgebra_glm::epsilon()
        && pb.abs() > nalgebra_glm::epsilon()
}

/// Finds the index of the vertex that is closest to the specified target vertex `t` among a set of
/// `points`.
pub fn find_closest_to(t: &Vec2, points: &Vec<Vec2>) -> (usize, f32) {
    // First, calculate the pair-wise distance between `t` and all of the other points
    let mut distances = points
        .iter()
        .map(|p| nalgebra_glm::distance(t, p))
        .collect::<Vec<_>>();

    // Then, find the index (and value) of the smallest such distance
    let mut min_distance = std::f32::MAX;
    let mut min_index = 0;

    for (index, d) in distances.iter().cloned().enumerate() {
        if d < min_distance {
            min_distance = d;
            min_index = index;
        }
    }

    (min_index, min_distance)

    // let min_index_and_val = distances
    //     .iter()
    //     .cloned()
    //     .enumerate()
    //     .min_by(|a, b| a.1.partial_cmp(&b.1).expect("Tried to compare a NaN"));
    //
    // min_index_and_val.unwrap()
}

/// Determines whether the specified points are "close to" each other, i.e. the same.
pub fn close_to(a: &Vec2, b: &Vec2) -> bool {
    nalgebra_glm::distance(a, b) <= nalgebra_glm::epsilon()
}

/// Calculates the (sub)set of "unique" points from an input set, within epsilon.
pub fn unique_points_among(ps: &Vec<Vec2>) -> Vec<Vec2> {
    let mut unique = vec![];

    for (i, a) in ps.iter().enumerate() {
        // Assume that this point `a` is unique
        let mut is_unique = true;

        // Iterate through all other points and check if `a` ~= `b` (taking care not
        // to compare `a` with itself in the inner for-loop)
        for (j, b) in ps.iter().enumerate() {
            if i != j && close_to(a, b) {
                is_unique = false;
                break;
            }
        }

        if is_unique {
            unique.push(*a);
        }
    }

    unique
}

/// Calculates the point of intersection between two infinite lines defined by the segments
/// `ab` and `cd`. The point of intersection might not lie within either of these segments
/// but elsewhere along their respective lines.
pub fn calculate_line_intersection(a: &Vec2, b: &Vec2, c: &Vec2, d: &Vec2) -> Option<Vec2> {
    let determinant = (d.y - c.y) * (b.x - a.x) - (d.x - c.x) * (b.y - a.y);

    if determinant.abs() > nalgebra_glm::epsilon() {
        let u_a = ((d.x - c.x) * (a.y - c.y) - (d.y - c.y) * (a.x - c.x)) / determinant;
        let u_b = ((b.x - a.x) * (a.y - c.y) - (b.y - a.y) * (a.x - c.x)) / determinant;

        if !((0.0 <= u_a && u_a <= 1.0) && (0.0 <= u_b && u_b <= 1.0)) {
            return None;
        }

        let x = a.x + u_a * (b.x - a.x);
        let y = a.y + u_a * (b.y - a.y);

        return Some(Vec2::new(x, y));
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple() {
        // See: https://brilliant.org/wiki/triangles-incenter/

        // Non-degenerate triangle
        let a = Vec2::new(0.0, 0.0);
        let b = Vec2::new(14.0, 0.0);
        let c = Vec2::new(5.0, 12.0);
        let res = calculate_triangle_incenter(&a, &b, &c);
        assert_eq!(Vec2::new(6.0, 4.0), res.unwrap().0);

        // Degenerate triangle
        let c = Vec2::new(18.0, 0.0); // Somewhere along x-axis, parallel to side AB
        let res = calculate_triangle_incenter(&a, &b, &c);
        assert_eq!(None, res);

        // Parallel lines
        let a = Vec2::new(0.0, 0.0);
        let b = Vec2::new(5.0, 0.0);
        let c = Vec2::new(0.0, 0.0);
        let d = Vec2::new(-5.0, 0.0);

        // Opposite directions (but parallel)
        assert_eq!(true, parallel(&a, &b, &c, &d));

        // Same directions (and parallel)
        assert_eq!(true, parallel(&a, &b, &d, &c));
    }
}

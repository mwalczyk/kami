use nalgebra_glm::{scale, Vec2};
use std::thread::current;

/// Returns the signed area of a triangle with vertices `a`, `b`, and `c`. In
/// particular, the area will be negative if the vertices are ordered in a CW
/// fashion and positive otherwise.
pub fn area2(a: &Vec2, b: &Vec2, c: &Vec2) -> f32 {
    (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)
}

/// Returns `true` if the point `c` is to the left of the oriented segment `ab`
/// and `false` otherwise.
pub fn left(a: &Vec2, b: &Vec2, c: &Vec2) -> bool {
    area2(a, b, c) > 0.0
}

/// Returns `true` if the point `c` is to the left of (or exactly on) the oriented
/// segment `ab` and `false` otherwise.
pub fn left_on(a: &Vec2, b: &Vec2, c: &Vec2) -> bool {
    area2(a, b, c) >= 0.0
}

/// Returns `true` if the point `c` is exactly on segment `ab` and `false` otherwise.
pub fn collinear(a: &Vec2, b: &Vec2, c: &Vec2) -> bool {
    area2(a, b, c) == 0.0
}

/// Determines whether two line segments `ab` and `cd` intersect. A "proper" intersection
/// does not include the case where the endpoint of one segment (say `c`) lies somewhere
/// on the other (closed) segment `ab`. For this, see `intersect()`.
pub fn intersect_prop(a: &Vec2, b: &Vec2, c: &Vec2, d: &Vec2) -> bool {
    // First, rule out the case where one (or both) endpoints of a segment is collinear
    // with the other
    if collinear(a, b, c) || collinear(a, b, d) || collinear(c, d, a) || collinear(c, d, b) {
        return false;
    }

    // Use XOR here to ensure that the points `c` and `d` are on opposite sides of
    // segment `ab` and vice-versa
    // println!("\t left(a, b, c): {}", left(a, b, c));
    // println!("\t left(a, b, d): {}", left(a, b, d));
    // println!("\t left(c, d, a): {}", left(c, d, a));
    // println!("\t left(c, d, b): {}", left(c, d, b));
    (left(a, b, c) ^ left(a, b, d)) && (left(c, d, a) ^ left(c, d, b))
}

/// Determines whether a point `c` is between points `a` and `b`. `c` is said to be between `a` and
/// `b` if it is collinear (i.e. on) the line containing segment `ab` but "between" `a` and `b` (not
/// beyond/before). Note that this function will return `true` if `c` is coincident to `a` or `b`.
pub fn between(a: &Vec2, b: &Vec2, c: &Vec2) -> bool {
    // `c` must be collinear to segment `ab` if it is between `a` and `b`
    if !collinear(a, b, c) {
        return false;
    }

    // If `ab` is not vertical, check "between-ness" based on the x-coordinate - otherwise,
    // check based on the y-coordinate
    if a.x != b.x {
        return ((a.x <= c.x) && (c.x <= b.x)) || ((a.x >= c.x) && (c.x >= b.x));
    }

    // `ab` is vertical
    ((a.y <= c.y) && (c.y <= b.y)) || ((a.y >= c.y) && (c.y >= b.y))
}

/// Determines whether the line segments `ab` and `cd` intersect. This function returns `true` if
/// the two segments intersect properly or if one endpoint of one segment lies between the two
/// endpoints of the other segment.
pub fn intersect(a: &Vec2, b: &Vec2, c: &Vec2, d: &Vec2) -> bool {
    if intersect_prop(a, b, c, d) {
        //println!("\tProper intersection found");
        return true;
    }

    // Is one of the endpoints of one segment between (i.e. contained in) the endpoints of the other?
    if between(a, b, c) || between(a, b, d) || between(c, d, a) || between(c, d, b) {
        //println!("\tIntersection found because one endpoint lies exactly on the other segment");
        return true;
    }

    false
}

#[derive(Debug)]
pub enum PolygonError {
    Degenerate,
}

#[derive(Clone, Debug)]
pub struct Polygon {
    vertices: Vec<Vec2>,
}

impl Polygon {
    /// Constructs a new polygon from the specified set of vertices.
    pub fn new(vertices: &Vec<Vec2>) -> Result<Polygon, PolygonError> {
        // Any valid polygon must have at least 3 vertices
        if vertices.len() < 3 {
            return Err(PolygonError::Degenerate);
        }

        let mut polygon = Polygon {
            vertices: vertices.clone(),
        };

        // Ensure that the polygon's vertices are always ordered in a CCW fashion
        // (if they are not, the signed area will be negative)
        if polygon.area_poly2() < 0.0 {
            polygon.vertices = polygon.vertices.into_iter().rev().collect();
            println!("Polygon vertices were CW - reversing to make CCW");
        }

        Ok(polygon)
    }

    /// Constructs a new regular polygon with the specified number of sides, circumradius, and center.
    pub fn regular(sides: usize, circumradius: f32, center: &Vec2) -> Polygon {
        let mut vertices = vec![];
        for i in 0..sides {
            vertices.push(Vec2::new(
                center.x
                    + circumradius * (2.0 * std::f32::consts::PI * (i as f32 / sides as f32)).cos(),
                center.y
                    + circumradius * (2.0 * std::f32::consts::PI * (i as f32 / sides as f32)).sin(),
            ));
        }
        Polygon::new(&vertices).unwrap()
    }

    /// Returns an immutable reference to the polygon's vertices.
    pub fn get_vertices(&self) -> &Vec<Vec2> {
        &self.vertices
    }

    /// Uniformly scales the entire polygon by an amount `s`.
    pub fn uniform_scale(&mut self, s: f32) {
        for vertex in self.vertices.iter_mut() {
            *vertex = *vertex * s;
        }
    }

    /// Translates the entire polygon by a vector `t`.
    pub fn translate(&mut self, t: &Vec2) {
        for vertex in self.vertices.iter_mut() {
            *vertex = *vertex + t;
        }
    }

    /// Returns the center of the polygon, i.e. the average of all of its vertices.
    pub fn center(&self) -> Vec2 {
        self.vertices.iter().sum::<Vec2>() / (self.vertices.len() as f32)
    }

    /// Computes twice the area of the polygon. Used internally to calculate the area of the
    /// polygon (see `area()`).
    fn area_poly2(&self) -> f32 {
        let mut sum = 0.0;

        // The "head" (fixed) vertex: arbitrarily chosen to be the first vertex
        let a = self.vertices[0];

        // Iterate over all other pairs of vertices in the polygon
        for i in 1..self.vertices.len() - 1 {
            let b = self.vertices[i + 0];
            let c = self.vertices[i + 1];
            sum += area2(&a, &b, &c);
        }

        sum
    }

    /// Computes the area of the polygon.
    pub fn area(&self) -> f32 {
        self.area_poly2().abs() * 0.5
    }

    /// Returns the two vertices that are neighbors to the vertex at `index`, wrapping around
    /// in a cyclical fashion, as necessary.
    fn get_neighbor_indices(&self, index: usize) -> (usize, usize) {
        let neighbor_prev = if index == 0 {
            self.vertices.len() - 1
        } else {
            index - 1
        };

        let neighbor_next = if index == self.vertices.len() - 1 {
            0
        } else {
            index + 1
        };

        (neighbor_prev, neighbor_next)
    }

    /// Determines whether or not the vertices at indices `i` and `j` are neighbors.
    fn are_neighbors(&self, i: usize, j: usize) -> bool {
        let (prev, next) = self.get_neighbor_indices(i);
        prev == j || next == j
    }

    /// Determines whether the vertex at index `other` is inside the cone defined by the
    /// vertex at index `apex`.
    fn in_cone(&self, apex: usize, other: usize) -> bool {
        let a = self.vertices[apex];
        let b = self.vertices[other];

        let (prev, next) = self.get_neighbor_indices(apex);
        let a0 = self.vertices[prev];
        let a1 = self.vertices[next];

        // First, check if `a` is a convex vertex
        if left_on(&a, &a1, &a0) {
            // Check if `a0` is strictly left of the segment `ab` and `a1` is
            // strictly right of the segment `ab`: if both of these conditions
            // are true, then `b` is inside the cone of vertex `a`
            return left(&a, &b, &a0) && left(&b, &a, &a1);
        }

        // At this point, we know that `a` is reflex: check if `a1` is left (or
        // on) segment `ab` and `a0` is right (or on) segment `ab`: if both of these
        // conditions are true, then `b` is inside the *inverted* cone of vertex `a`
        // (i.e. the cone if we consider the complement of reflex vertex `a`)
        //
        // Therefore, we negate this result to check if `b` is in the interior
        !(left_on(&a, &b, &a1) && left_on(&b, &a, &a0))
    }

    /// Checks whether the line segment connecting the vertices at indices `i` and `j`
    /// form a diagonal, i.e. a segment that does not intersect any of the other edges
    /// of the polygon and lies completely in the interior of the polygon.
    pub fn diagonal(&self, i: usize, j: usize) -> bool {
        // A single vertex can't form a diagonal
        if i == j {
            return false;
        }

        // See if the diagonal formed between the i-th and j-th vertices of this polygon
        // intersects any edges that don't contain either of the diagonal's endpoints
        for k in 0..self.vertices.len() {
            if k != i
                && (k + 1) % self.vertices.len() != i
                && k != j
                && (k + 1) % self.vertices.len() != j
                && intersect(
                    &self.vertices[i],
                    &self.vertices[j],
                    &self.vertices[k],
                    &self.vertices[(k + 1) % self.vertices.len()],
                )
            {
                return false;
            }
        }

        // At this point, we know that the diagonal doesn't intersect any of the polygon
        // edges, but we still need to check that it lies entirely in the interior of the
        // polygon
        self.in_cone(i, j) && self.in_cone(j, i)
    }

    pub fn triangulate(&self) -> Vec<Polygon> {
        // This polygon is a triangle, so we are done
        if self.vertices.len() == 3 {
            return vec![self.clone()];
        }

        let mut copy = self.clone();

        // First, iterate through each vertex and determine whether it is an "ear": a
        // vertex is an "ear" if the segment connecting the vertex's prev/next neighbors
        // is a valid diagonal of this polygon
        let mut is_ear = vec![];

        for i in 0..self.vertices.len() {
            let (prev, next) = copy.get_neighbor_indices(i);
            is_ear.push(copy.diagonal(prev, next));
        }
        assert_eq!(copy.vertices.len(), is_ear.len());
        println!("{:?}", is_ear);

        let mut triangles = vec![];

        while copy.vertices.len() > 3 {

            for v2 in 0..copy.vertices.len() {

                if is_ear[v2] {
                    // Find 5 consecutive vertices, with v2 at the center: in CCW order, we have
                    // v0, v1, v2, v3, v4
                    let (mut v1, mut v3) = copy.get_neighbor_indices(v2);
                    let (_, mut v4) = copy.get_neighbor_indices(v3);
                    let (mut v0, _) = copy.get_neighbor_indices(v1);

                    // This ear forms a triangle
                    triangles.push(
                        Polygon::new(&vec![
                            copy.vertices[v1],
                            copy.vertices[v2],
                            copy.vertices[v3],
                        ])
                        .unwrap(),
                    );

                    // Remove the triangle (i - 1, i, i + 1), which is an ear
                    copy.vertices.remove(v2);
                    is_ear.remove(v2);

                    // The indices that come after the removed vertex v2 will have been shifted
                    // over one place because of the removal operation: note that vertices v0 and v1
                    // may also need to be shifted
                    //
                    // Consider a square polygon where v2 is the first vertex. Then v1 will be the 3rd
                    // vertex and v0 will be the 2nd. So, we need to subtract 1 from these indices as
                    // well.
                    if v0 > v2 {
                        v0 -= 1;
                    }
                    if v1 > v2 {
                        v1 -= 1;
                    }
                    if v3 > v2 {
                        v3 -= 1;
                    }
                    if v4 > v2 {
                        v4 -= 1;
                    }

                    // Update ear status of the two vertices that were part of the ear that we
                    // just removed
                    is_ear[v1] = copy.diagonal(v0, v3);
                    is_ear[v3] = copy.diagonal(v1, v4);

                    break;
                }
            }
        }

        // At this point, the remaining polygon will necessarily be a triangle, so we
        // need to make sure to add it as well
        triangles.push(copy);

        triangles
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple() {
        let polygon = Polygon::new(&vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(5.0, 0.0),
            Vec2::new(0.0, 5.0),
        ])
        .unwrap();

        let area = polygon.area();
        println!("Polygon area is: {:?}", area);
        approx::assert_ulps_eq!(12.5, area);

        println!("Testing left, collinear, between, etc.");

        // Vertical
        let a = Vec2::new(0.0, 0.0);
        let b = Vec2::new(0.0, 1.0);
        let c = Vec2::new(0.0, 0.5);
        let d = Vec2::new(-1.0, 0.0);
        let e = Vec2::new(1.0, 0.0);
        assert!(between(&a, &b, &c));
        assert!(left(&a, &b, &d));
        assert!(!left(&a, &b, &e));

        // Non-vertical
        let a = Vec2::new(0.0, 0.0);
        let b = Vec2::new(1.0, 1.0);
        let c = Vec2::new(0.5, 0.5);
        let d = Vec2::new(0.0, 1.0);
        let e = Vec2::new(1.0, 0.0);
        assert!(between(&a, &b, &c));
        assert!(left(&a, &b, &d));
        assert!(!left(&a, &b, &e));

        // Square
        // let polygon = Polygon::new(&vec![
        //     Vec2::new(0.0, 0.0),
        //     Vec2::new(1.0, 0.0),
        //     Vec2::new(1.0, 1.0),
        //     Vec2::new(0.0, 1.0),
        // ]).unwrap();
        //
        // let triangles = polygon.triangulate();
        // println!("{:?}", triangles);

        // From book
        let polygon = Polygon::new(&vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(10.0, 7.0),
            Vec2::new(12.0, 3.0),
            Vec2::new(20.0, 8.0),
            Vec2::new(13.0, 17.0),
            Vec2::new(10.0, 12.0),
            Vec2::new(12.0, 14.0),
            Vec2::new(14.0, 9.0),
            Vec2::new(8.0, 10.0),
            Vec2::new(6.0, 14.0),
            Vec2::new(10.0, 15.0),
            Vec2::new(7.0, 18.0),
            Vec2::new(0.0, 16.0),
            Vec2::new(1.0, 13.0),
            Vec2::new(3.0, 15.0),
            Vec2::new(5.0, 8.0),
            Vec2::new(-2.0, 9.0),
            Vec2::new(5.0, 5.0),
        ])
        .unwrap();

        let mut is_ear = vec![];
        for i in 0..polygon.vertices.len() {
            let (prev, next) = polygon.get_neighbor_indices(i);
            //println!("Examining neighbors of {}, which are {} and {}", i, prev, next);

            is_ear.push(polygon.diagonal(prev, next));
        }
        println!("{:?}", is_ear);

        // println!("is 14 left of 11->12? {}", left(&polygon.vertices[11], &polygon.vertices[12], &polygon.vertices[14]));
        // println!("is 16 left of 11->12? {}", left(&polygon.vertices[11], &polygon.vertices[12], &polygon.vertices[16]));

        //  println!("does 14->16 intersect 11->12? {}", intersect(&polygon.vertices[14], &polygon.vertices[16], &polygon.vertices[11], &polygon.vertices[12]));

        let triangles = polygon.triangulate();
        println!("Triangulation resulted in {:?} triangles", triangles.len());
    }
}

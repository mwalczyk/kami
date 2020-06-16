pub enum SectorType {
    Anto,
    Iso,
}

pub enum MajorityType {
    MountainLike,
    ValleyLike,
}

struct Vertex {}

impl Vertex {
    /// The degree of a vertex is the number of creases incident.
    fn degree() -> usize {
        unimplemented!()
    }

    /// A vertex is developable if the sum of its sector angles is 360.
    fn is_developable() -> bool {
        unimplemented!()
    }

    /// Returns the sector angles surrounding this vertex, sorted in CCW order.
    fn sector_angles() -> Vec<f32> {
        unimplemented!()
    }

    /// Returns this vertex's majority type (i.e. whether it is surrounded by
    /// more M creases or V creases). If the vertex is not flat-foldable
    /// (because it does not satisfy the Maekawa-Justin Theorem or one or more
    /// of its incident creases are unassigned), `None` will be returned.
    fn majority_type() -> Option<MajorityType> {
        unimplemented!()
    }

    /// The vertex type is a string of characters composed of Ms and Vs giving
    /// the fold types that one encounters as one goes around the vertex CCW. For
    /// example, `MVMVVM`.
    fn vertex_type() -> String {
        unimplemented!()
    }

    /// Returns `true` if this is an interior vertex (i.e. located somewhere
    /// in the interior of the paper) and `false` otherwise.
    fn is_interior() -> bool {
        unimplemented!()
    }
}

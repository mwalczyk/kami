pub enum Assignment {
    Mountain,
    Valley,
    Border,
    Unfolded,
}

pub struct Crease {}

impl Crease {
    fn dihedral_angle() -> f32 {
        unimplemented!()
    }

    fn fold_angle() -> f32 {
        unimplemented!()
    }

    /// The angle that this crease makes with respect to some reference line
    /// (for example, the positive `x`-axis).
    fn fold_direction_angle() -> f32 {
        unimplemented!()
    }
}

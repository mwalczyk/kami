pub trait InnerSpace {
    fn dot();

    fn normalize();
}

pub trait OuterSpace {
    fn wedge();
}

//https://github.com/rustgd/cgmath/blob/master/src/structure.rs
//https://crates.parity.io/alga/linear/trait.InnerSpace.html

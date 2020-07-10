[![Build Status](https://travis-ci.org/mwalczyk/kami.svg?branch=master)](https://travis-ci.org/github/mwalczyk/kami) 

# Kami

🦢 A library for computational origami and geometric algebra (vectors, bivectors, rotors, etc.).

## Resources

### Geometric Algebra
- [Clifford / Geometric Algebra Intuition](https://www.euclideanspace.com/maths/algebra/clifford/index.htm)
- [Geometric Algebra: An Introduction with Applications in Euclidean and Conformal Geometry](https://scholarworks.sjsu.edu/cgi/viewcontent.cgi?article=7943&context=etd_theses)
- [Geometric Algebra Primer](http://www.jaapsuter.com/geometric-algebra.pdf)
- [Bivector (and Ganja.js)](https://bivector.net/doc.html)
- [Exterior Algebra](https://en.wikipedia.org/wiki/Exterior_algebra)
- [Visualizing Geometric Product Relationships](https://www.shapeoperator.com/2019/01/07/relating-dot-wedge/)
- [An Introduction to Geometric Algebra](https://bitworking.org/news/ga/2d/)

### YouTube Series
See the following [playlist](https://www.youtube.com/playlist?list=PLpzmRsG7u_gqaTo_vEseQ7U8KFvtiJY4K).

### Other Implementations
- [Marc Ten Bosch](https://marctenbosch.com/quaternions/code.htm)
- [Versor C++11 Library from UCSB](http://versor.mat.ucsb.edu/)
- [Ultraviolet Rust Library](https://github.com/termhn/ultraviolet)

### Misc
- [Generics vs. Associated Types in Rust](https://stackoverflow.com/questions/32059370/when-is-it-appropriate-to-use-an-associated-type-versus-a-generic-type)

### WASM
- [Getting Started](https://dev.to/sendilkumarn/rust-and-webassembly-for-the-masses-wasm-pack-3d6p)
- [Yew](https://yew.rs/docs/)


### Graph Data Structures in Rust 
- [Graphlib](https://github.com/purpleprotocol/graphlib)
- [Modeling Graphs in Rust Using Indices](http://smallcultfollowing.com/babysteps/blog/2015/04/06/modeling-graphs-in-rust-using-vector-indices/)
- [Petgraph](https://github.com/petgraph/petgraph)

## Roadmap

- Math
	- Add remaining tests to all modules
	- Potentially remove the trivector module
	- Sort out inner / outer / norm traits
	- Verify correctness of all algorithms
	- Provide some sort of documentation / geometric interpretations
- Graphs
	- Planar graphs (see [petgraph](https://docs.rs/petgraph/0.5.1/petgraph/))
	- Graphviz output format
	- .FOLD import / export
- Origami
	- Implement HH axioms using geometric algebra
- Geometry
	- Read "Computational Geometry in C" and implement polygon class
	- (Re)implement robust intersection routines
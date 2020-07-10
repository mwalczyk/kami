[![Build Status](https://travis-ci.org/mwalczyk/kami.svg?branch=master)](https://travis-ci.org/github/mwalczyk/kami) 

# Kami

🦢 A (very) W.I.P. library for computational origami.

## Description

After studying computational origami for the last 6-8 months and making a variety of rather disjoint projects, I am working towards creating a unified codebase for future origami work. This is *very* much a work-in-progress, and it will be a few months before the library is even remotely usable. It's my hope that this repo becomes a sort of "diary" / "scratch pad" for future studies: I will add more robustness and functionality as I continue to learn more about origami design, kinematics, and computational geometry. Pretty much all of the ideas here are taken from Robert Lang's books, _Twists, Tilings, and Tessellations_ and/or _Origami Design Secrets_, so I claim no ownership over any of this knowledge. 

Initially, I set out to build this library using geometric (or Clifford) algebra: a unified mathematical framework that includes "traditional" linear algebra, quaternions, complex numbers, and many other mathematical phenomena, all in one "package." This is still a goal of mine, but for now, I've decided to focus on some of the other higher-level functionality, deferring to the already wonderful `nalgebra` library for basic mathematical operations.

Because the Rust GUI story isn't great yet (although several libraries, like `Iced` and `Druid` are getting there), I've opted to create visualizations using the browser. This pretty straightforward with WASM and the `websys` crate. Regardless, this library is meant to be used _within_ other projects that have a graphical interface - it isn't a front-end application itself. 

In other words, you would use this library internally to build things like:
- A crease pattern editor
- A tessellation generator
- Algorithmic origami designs
- Rigid origami simulations

I'm using Rust because of its awesome type / ownership system, tooling, and community 🦀. 

## References

- Geometric Algebra
  - [Clifford / Geometric Algebra Intuition](https://www.euclideanspace.com/maths/algebra/clifford/index.htm)
  - [Geometric Algebra: An Introduction with Applications in Euclidean and Conformal Geometry](https://scholarworks.sjsu.edu/cgi/viewcontent.cgi?article=7943&context=etd_theses)
  - [Geometric Algebra Primer](http://www.jaapsuter.com/geometric-algebra.pdf)
  - [Bivector (and Ganja.js)](https://bivector.net/doc.html)
  - [Exterior Algebra](https://en.wikipedia.org/wiki/Exterior_algebra)
  - [Visualizing Geometric Product Relationships](https://www.shapeoperator.com/2019/01/07/relating-dot-wedge/)
  - [An Introduction to Geometric Algebra](https://bitworking.org/news/ga/2d/)
  - See the following [playlist](https://www.youtube.com/playlist?list=PLpzmRsG7u_gqaTo_vEseQ7U8KFvtiJY4K).
  - [Marc Ten Bosch](https://marctenbosch.com/quaternions/code.htm)
  - [Versor C++11 Library from UCSB](http://versor.mat.ucsb.edu/)
  - [Ultraviolet Rust Library](https://github.com/termhn/ultraviolet)
- Rust / WASM
  - [Generics vs. Associated Types in Rust](https://stackoverflow.com/questions/32059370/when-is-it-appropriate-to-use-an-associated-type-versus-a-generic-type)
  - [Getting Started](https://dev.to/sendilkumarn/rust-and-webassembly-for-the-masses-wasm-pack-3d6p)
  - [Yew](https://yew.rs/docs/)
- Graphs
  - [Graphlib](https://github.com/purpleprotocol/graphlib)
  - [Modeling Graphs in Rust Using Indices](http://smallcultfollowing.com/babysteps/blog/2015/04/06/modeling-graphs-in-rust-using-vector-indices/)
  - [Petgraph](https://github.com/petgraph/petgraph)

## Roadmap

- Math
  - Vectors (2D and 3D)
  - Matrices (primarily 2x3, 3x3, 3x4, and 4x4)
  - Linear algebra operations
    - Moore-Penrose pseudo-inverse
    - Inverse
    - Transpose
    - Determinant
  - Revisit geometric / Clifford algebra 
  - Implement HH axioms using geometric algebra
- Graphs
  - Planar graphs
    - Add new vertex (while maintaining planarity)
    - Add new edge (while maintaining planarity)
    - Delete an existing vertex (and all incident edges - optional)
    - Delete an existing edge (and any stray vertices - optional)
    - Determine 2-colorability
    - Graph neural networks (GNNs)?
    - Search (for a particular node, edge, or face)
    - Traversal
      - Grab all reference points (i.e. nodes)
      - Grab all fold vectors (i.e. vectors oriented along each edge of the graph)
      - Grab all interior fold intersections with outward facing, sorted edges emanating from each
      - Grab face boundaries (i.e. CCW sorted edges that form each face)
      - Angle of each fold vector w.r.t. +x-axis
      - Grab all face corner angles
      - Grab all face corner points
      - Construct fold paths
      - Construct fold map (matrix transformations accumulated by traversing a fold path)
  - I/O formats
    - Graphviz output format (.dot)
    - .FOLD import / export
    - .oripa import / export
    - .svg import / export
- Origami
  - Crease
    - Specify assignment
    - Specify fold angle (target, upper bound, lower bound)
  - Crease pattern
    - Maintains a (mutable) planar graph
    - Compute the stacking order (Oripa "translucent" preview)
    - Calculate vertex properties
      - Degree
      - Iso / anto
      - Majority type (mountain-like or valley-like)
      - Sector angles
    - Flat-foldability checks
      - Big-little-big-angle theorem (BLBA)
      - Kawasaki-Justin theorem (KJT)
      - Algebraic KJT (AKJT) - same as KJT but with matrices
      - Maekawa-Justin theorem (MJT)
      - Single vertex flat-foldable test (SVFFT) - accumulation of the previous 3 theorems
    - Calculate Justin path (similar to fold map from "Active Origami"?)
    - Apply a sequence of Huzita-Hatori axioms, starting from a blank sheet of "paper"
  - Misc.
    - Universal molecule
    - Tree theory
    - Box pleating
    - 22.5 design
- Geometry
  - Read "Computational Geometry in C" 
  - Polygons (regular, irregular)
    - Compute area
    - Compute perimeter
    - Compute interior angles
    - Compute exterior angles
    - Apply transforms (rotate, scale, translate, etc.)
    - Calculate incenter, circumradius, etc.
    - Calculate bounding box
    - Degeneracy tests
    - Calculate angle bisectors
    - "Slicing" operations (fracture polygon into one or more sub-polygons)
    - Polygon to polyline conversion
    - Calculate "skeleton" (see ODS)
  - Lines (segments, infinite)
  - Polylines
    - Add / remove points
    - Resample / refine operations
    - Smoothing operations
    - Calculate length
    - Calculate tangents
  - Tilings (Archimedean)
    - Vertex figures
    - Lattice vectors
    - Array of polygons with offsets and rotations
  - Algorithms / utilities
    - Triangulate an arbitrary polygon
    - Compute the convex hull of a set of 2D points
    - Circle packing (bindings to Python SciPy (SLSQP) optimization)
    - Intersection routines
      - Line-line
      - Line-polygon
      - Segment-segment
      - Segment-polygon
      - Polygon-polygon
    - Misc.
      - Perpendicular bisector
      - Point-on-segment test
      - Point-on-line test
      - Point adjacency 
      - Point sorting
      
    

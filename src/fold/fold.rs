use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

/// A .FOLD file format parser. Note that multiple frames are not yet supported.
/// See: https://github.com/edemaine/fold
#[derive(Serialize, Deserialize, Debug)]
pub struct Fold {
    /// File metadata

    /// The version of the FOLD spec that the file assumes (a number)
    #[serde(default)]
    pub file_spec: i32,

    /// The software that created the file (a string)
    #[serde(default)]
    pub file_creator: String,

    /// The human author (a string)
    #[serde(default)]
    pub file_author: String,

    /// A title for the entire file (a string)
    #[serde(default)]
    pub file_title: String,

    /// A description of the entire file (a string)
    #[serde(default)]
    pub file_description: String,

    /// A subjective interpretation about what the entire file represents (array of strings)
    #[serde(default)]
    pub file_classes: Vec<String>,

    #[serde(skip)]
    pub file_frames: String,

    /// Frame metadata

    /// The human author (a string)
    #[serde(default)]
    pub frame_author: String,

    /// A title for the frame (a string)
    #[serde(default)]
    pub frame_title: String,

    /// A description of the frame (a string)
    #[serde(default)]
    pub frame_description: String,

    /// A subjective interpretation about what the frame represents (array of strings)
    #[serde(default)]
    pub frame_classes: Vec<String>,

    /// Attributes that objectively describe properties of the folded structure being represented (array of strings)
    #[serde(default)]
    pub frame_attributes: Vec<String>,

    /// Physical or logical unit that all coordinates are relative to (a string)
    #[serde(default)]
    pub frame_unit: String,

    /// Vertex information

    /// For each vertex, an array of coordinates, such as [x, y, z] or [x, y] (where z is implicitly zero)
    #[serde(default)]
    pub vertices_coords: Vec<Vec<f32>>,

    /// For each vertex, an array of vertices (vertex IDs) that are adjacent along edges
    #[serde(default)]
    pub vertices_vertices: Vec<Vec<i32>>,

    /// For each vertex, an array of face IDs for the faces incident to the vertex
    #[serde(default)]
    pub vertices_faces: Vec<Vec<i32>>,

    /// Edge information

    /// For each edge, an array [u, v] of two vertex IDs for the two endpoints of the edge
    #[serde(default)]
    pub edges_vertices: Vec<Vec<i32>>,

    /// For each edge, an array of face IDs for the faces incident to the edge
    #[serde(default)]
    pub edges_faces: Vec<Vec<i32>>,

    /// For each edge, a string representing its fold direction assignment
    #[serde(default)]
    pub edges_assignment: Vec<String>,

    /// For each edge, the fold angle (deviation from flatness) along each edge of the pattern
    #[serde(default, rename(deserialize = "edges_foldAngles"))]
    pub edges_foldangles: Vec<f32>,

    /// For each edge, the length of the edge (a number)
    #[serde(default)]
    pub edges_lengths: Vec<f32>,

    /// Face information

    /// For each face, an array of vertex IDs for the vertices around the face in counterclockwise order
    #[serde(default)]
    pub faces_vertices: Vec<Vec<i32>>,

    /// For each face, an array of edge IDs for the edges around the face in counterclockwise order
    #[serde(default)]
    pub faces_edges: Vec<Vec<i32>>,

    /// Layer information

    /// An array of triples [f, g, s] where f and g are face IDs and s is an integer between −1 and 1
    #[serde(default, rename(deserialize = "faceOrders"))]
    pub faceorders: Vec<[i32; 3]>,

    /// An array of triples [e, f, s] where e and f are edge IDs and s is an integer between −1 and 1
    #[serde(default, rename(deserialize = "edgeOrders"))]
    pub edgeorders: Vec<[i32; 3]>,
}

impl Fold {
    //pub fn new() -> Fold {}

    pub fn from_file(path: &Path) -> std::io::Result<Fold> {
        let mut file = File::open(path)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;

        let mut fold: Fold = serde_json::from_str(&contents).expect("Failed to load JSON");
        fold.resolve();
        fold.validate();

        Ok(fold)
    }

    pub fn resolve(&mut self) {
        // Fill in missing bits of information from .fold spec
    }

    pub fn validate(&self) -> bool {
        // vertices_coords, vertices_vertices, and vertices_faces should all be the same length (or zero), etc.
        true
    }
}

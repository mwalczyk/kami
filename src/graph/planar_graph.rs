use nalgebra_glm::Vec2;

pub type NodeIndex = usize;
pub type EdgeIndex = usize;
pub type FaceIndex = usize;

const DEFAULT_LABEL: &str = "";

pub struct Node<T> {
    /// The 2D position of this node
    position: Vec2,

    /// Data associated with this node
    data: T,

    /// Node label, useful for debugging purposes
    label: String,
}

impl<T> Node<T> {
    pub fn new(position: &Vec2, data: T, label: &str) -> Node<T> {
        Node {
            position: *position,
            data,
            label: String::from(label),
        }
    }

    pub fn get_position(&self) -> &Vec2 {
        &self.position
    }
}

pub struct Edge<T> {
    /// The index of the first vertex incident to this edge
    src: NodeIndex,

    /// The index of the second vertex incident to this edge
    dst: NodeIndex,

    /// Data associated with this edge (crease assignment, for example)
    data: T,

    /// Edge label, useful for debugging purposes
    label: String,
}

impl<T> Edge<T> {
    pub fn new(src: NodeIndex, dst: NodeIndex, data: T, label: &str) -> Edge<T> {
        Edge {
            src,
            dst,
            data,
            label: String::from(label),
        }
    }

    pub fn get_indices(&self) -> (NodeIndex, NodeIndex) {
        (self.src, self.dst)
    }
}

pub enum PlanarGraphError {
    NoSuchVertex,
    NoSuchEdge,
}

/// A planar graph data structure with generic node data `T` and edge data `U`.
pub struct PlanarGraph<T, U> {
    nodes: Vec<Node<T>>,
    edges: Vec<Edge<U>>,
}

impl<T, U> PlanarGraph<T, U> {
    pub fn new() -> Self {
        PlanarGraph {
            nodes: vec![],
            edges: vec![],
        }
    }

    pub fn node_count(self) -> usize {
        self.nodes.len()
    }

    pub fn edge_count(self) -> usize {
        self.edges.len()
    }

    pub fn get_nodes(&self) -> &Vec<Node<T>> {
        &self.nodes
    }

    pub fn get_edges(&self) -> &Vec<Edge<U>> {
        &self.edges
    }

    pub fn get_vertex(&self, id: NodeIndex) -> Option<&Node<T>> {
        self.nodes.get(id)
    }

    pub fn get_edge(&self, id: EdgeIndex) -> Option<&Edge<U>> {
        self.edges.get(id)
    }

    pub fn has_vertex(&self, id: NodeIndex) -> bool {
        id > 0 && id < self.nodes.len()
    }

    pub fn has_edge(&self, a: NodeIndex, b: NodeIndex) -> bool {
        for edge in self.edges.iter() {
            // Order does not matter in this graph implementation (it is undirected)
            if (a, b) == (edge.src, edge.dst) || (b, a) == (edge.src, edge.dst) {
                return true;
            }
        }
        false
    }

    /// Returns the number of edges that touch (i.e. include) the specified node
    pub fn number_of_edges_incident_to(&self, id: NodeIndex) -> Option<usize> {
        if let Some(node) = self.get_vertex(id) {
            return Some(
                self.edges
                    .iter()
                    .filter(|e| e.src == id || e.dst == id)
                    .count(),
            );
        }
        None
    }

    pub fn add_node(&mut self, position: &Vec2, data: T) {
        self.nodes.push(Node::new(position, data, DEFAULT_LABEL));
    }

    pub fn add_edge(
        &mut self,
        src: NodeIndex,
        dst: NodeIndex,
        data: U,
    ) -> Result<(), PlanarGraphError> {
        // Don't add duplicate edges
        if !self.has_edge(src, dst) {
            // Make sure `src` and `dst` vertices actually exist
            if let (Some(_), Some(_)) = (self.get_vertex(src), self.get_vertex(dst)) {
                self.edges.push(Edge::new(src, dst, data, DEFAULT_LABEL));
                return Ok(());
            }
            return Err(PlanarGraphError::NoSuchVertex);
        }

        // The edge already exists in the graph, which isn't an error, but we don't
        // need to do anything else here
        Ok(())
    }
}

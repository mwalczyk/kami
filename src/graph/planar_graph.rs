use crate::geometry::utils;

use nalgebra_glm::Vec2;

pub type NodeIndex = usize;
pub type EdgeIndex = usize;
pub type FaceIndex = usize;

const DEFAULT_LABEL: &str = "";

#[derive(Clone, Debug)]
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

#[derive(Clone, Debug)]
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
#[derive(Clone)]
pub struct PlanarGraph<T, U> {
    nodes: Vec<Node<T>>,
    edges: Vec<Edge<U>>,
}

impl<T, U> PlanarGraph<T, U>
where
    T: Copy + Clone + Default,
    U: Copy + Clone + Default,
{
    pub fn new() -> Self {
        PlanarGraph {
            nodes: vec![],
            edges: vec![],
        }
    }

    /// Returns the number of nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the number of edges in the graph.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Returns an immutable reference to the nodes in the graph.
    pub fn get_nodes(&self) -> &Vec<Node<T>> {
        &self.nodes
    }

    /// Returns an immutable reference to the edges in the graph.
    pub fn get_edges(&self) -> &Vec<Edge<U>> {
        &self.edges
    }

    /// Returns a reference to the node at index `id` if it exists and `None` otherwise.
    pub fn get_vertex(&self, id: NodeIndex) -> Option<&Node<T>> {
        self.nodes.get(id)
    }

    /// Returns a reference to the edge at index `id` if it exists and `None` otherwise.
    pub fn get_edge(&self, id: EdgeIndex) -> Option<&Edge<U>> {
        self.edges.get(id)
    }

    pub fn has_node_with_index(&self, id: NodeIndex) -> bool {
        id > 0 && id < self.nodes.len()
    }

    pub fn has_edge_with_index(&self, id: EdgeIndex) -> bool {
        id > 0 && id < self.edges.len()
    }

    pub fn has_edge_between(&self, a: NodeIndex, b: NodeIndex) -> bool {
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

    /// Attempts to add a new node to the graph at `position`. This function returns a tuple with 2 entries:
    ///
    /// 1. The index of the newly added node (or the index of an existing node if the specified node was
    ///    close to an existing node)
    /// 2. A list containing the indices of all of the edges that changed as a result of adding the new node
    pub fn add_node(&mut self, position: &Vec2, data: T) -> (NodeIndex, Vec<EdgeIndex>) {
        // Check if the vertex is the same as an existing vertex (within epsilon) - map the vec of node structs
        // into a vec of 2D positions so that we can call the appropriate function
        let node_positions = self
            .nodes
            .iter()
            .map(|node| node.position)
            .collect::<Vec<_>>();
        let (index, distance) = utils::find_closest_to(position, &node_positions);

        if distance < nalgebra_glm::epsilon() {
            // If the found distance is less than the specified threshold, the specified
            // node is considered a "duplicate," so we return the index of the existing
            // node
            println!("Attempting to add node that is very close to the node at index {} - returning the existing index instead ", index);
            return (index, vec![]);
        }

        // The new node is added at the end of the list, so we return that index along
        // with any edges that may have changed / been created
        self.nodes.push(Node::new(position, data, DEFAULT_LABEL));

        // Adding a new node along an existing edge(s) requires the offending edges to be split
        let changed_edges = self.split_edges_at_node(self.node_count() - 1);

        (self.node_count() - 1, changed_edges)
    }

    /// Attempts to add a new edge between the two specified node indices, splitting any intersecting
    /// edges and adding new nodes and edges as necessary to maintain planarity
    pub fn add_edge(
        &mut self,
        src: NodeIndex,
        dst: NodeIndex,
        data: U,
    ) -> Result<(Vec<NodeIndex>, Vec<EdgeIndex>), PlanarGraphError> {
        // Don't add duplicate (collinear) edges
        if !self.has_edge_between(src, dst) {
            // Make sure `src` and `dst` vertices actually exist
            if let (Some(_), Some(_)) = (self.get_vertex(src), self.get_vertex(dst)) {
                // First, push back the new edge
                self.edges.push(Edge::new(src, dst, data, DEFAULT_LABEL));

                // Perform line-segment/line-segment intersection tests
                let (changed_nodes, mut changed_edges) =
                    self.split_edges_along_edge(self.edge_count() - 1);

                // If the array of changed edges is empty, this means that no edge splitting was necessary,
                // but we still want to make sure we return at least one edge index (in this case, the edge
                // was simply added at the end of the graph's edge array, so just return that index)
                if changed_edges.is_empty() {
                    changed_edges.push(self.edge_count() - 1);
                }

                return Ok((changed_nodes, changed_edges));
            }

            return Err(PlanarGraphError::NoSuchVertex);
        }

        // The edge already exists in the graph, which isn't an error, but we don't
        // need to do anything else here
        Ok((vec![], vec![]))
    }

    /// Splits any edges that contain the specified node in their interiors. This function returns the
    /// indices of all newly split edges.
    fn split_edges_at_node(&mut self, id: NodeIndex) -> Vec<EdgeIndex> {
        println!("Splitting edges at node with index {}", id);

        let mut changed_edges = vec![];

        // TODO: check if the vertex at index `id` actually exists

        let mut additional_edges = vec![];

        for (index, edge) in self.edges.iter_mut().enumerate() {
            // Does the target node lie along this edge somewhere?
            let on_edge = utils::is_on_line_segment(
                &self.nodes[edge.get_indices().0].position,
                &self.nodes[edge.get_indices().1].position,
                &self.nodes[id].position,
                false,
            );

            if on_edge {
                println!("Node splits edge with index {}", index);

                // This edge isn't added to the graph's list of edges right away for 2 reasons:
                // 1. Rust's borrow checker
                // 2. We don't need to (re)test if the target node is on this edge, since we
                //    just created it
                additional_edges.push(Edge::new(
                    edge.get_indices().1,
                    id,
                    U::default(),
                    DEFAULT_LABEL,
                ));

                // Put one of the new edges in the slot that was previously occupied by the old,
                // un-split edge and the other at the end of the array - the prior is done so that
                // we don't have to worry about the rest of the edges being "shuffled" as a result
                // of a standard array "remove" operation...essentially, we want to add/remove edges
                // in-place
                *edge = Edge::new(edge.get_indices().0, id, U::default(), DEFAULT_LABEL);
                changed_edges.push(index);
            }
        }

        // Now, add the other edges that we constructed during splitting
        for edge in additional_edges.iter().cloned() {
            self.edges.push(edge);

            // The edge that we just added is at position `self.edges.len() - 1`, so add that index
            // to the list of changed edges
            changed_edges.push(self.edge_count() - 1);
        }

        changed_edges
    }

    fn split_edges_along_edge(&mut self, id: EdgeIndex) -> (Vec<NodeIndex>, Vec<EdgeIndex>) {
        // TODO: check if the edge at index `id` exists

        let mut changed_nodes = vec![];
        let mut changed_edges = vec![];

        // An array containing all of the points of intersections
        let mut intersections = vec![];

        for (index, edge) in self.edges.iter().enumerate() {
            // Don't intersect the target edge with itself
            if id != index {
                // The point of intersection (or null if no intersection is found)
                if let Some(intersection) = utils::calculate_line_segment_intersection(
                    &self.nodes[edge.get_indices().0].position,
                    &self.nodes[edge.get_indices().1].position,
                    &self.nodes[self.edges[id].get_indices().0].position,
                    &self.nodes[self.edges[id].get_indices().1].position,
                ) {
                    intersections.push(intersection);
                    println!("Found intersection between edge {} and edge {}", id, index);
                }
            }
        }

        // Get rid of duplicate intersection points - this happens, for example, when the target
        // edge passes through a node that is the shared endpoint of multiple neighbor edges:
        //
        //      \           <- existing edge A
        //       \
        //  ------*-----    <- target edge
        //       /
        //      /           <- existing edge B
        //
        // For example, in the diagram above, the pre-existing edges A and B terminate at a
        // common node, and when the target edge is added, TWO intersections will be found (one
        // against edge A, one against edge B), but the two points will coincide with one
        // another
        intersections = utils::unique_points_among(&intersections);
        println!(
            "Found {} unique intersections: {:?}",
            intersections.len(),
            intersections
        );

        for intersection in intersections.iter() {
            // One annoying edge-case we have to consider here is, what happens when the point of
            // intersection coincides (i.e. overlaps with) an existing node? In this case, we want
            // to split any offending edges along the existing node and continue - otherwise, the
            // point of intersection should be considered a new node, so we add it and split along
            // any edges, as necessary
            let node_positions = self
                .nodes
                .iter()
                .map(|node| node.position)
                .collect::<Vec<_>>();
            let (index, distance) = utils::find_closest_to(intersection, &node_positions);

            if distance < nalgebra_glm::epsilon() {
                println!("Point of intersection coincides with an existing node");
                changed_edges.extend(&self.split_edges_at_node(index));
            } else {
                // Add a new node at the point of intersection, which returns a node index and a list
                // of all of the edge indices that changed as a result of the additional node
                let (node_index, edge_indices) = self.add_node(intersection, T::default());

                changed_nodes.push(node_index);
                changed_edges.extend(&edge_indices);
            }
        }

        (changed_nodes, changed_edges)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple() {
        let mut graph: PlanarGraph<(), ()> = PlanarGraph::new();
        graph.add_node(&Vec2::new(100.0, 100.0), ());
        graph.add_node(&Vec2::new(200.0, 100.0), ());
        graph.add_node(&Vec2::new(150.0, 50.0), ());
        graph.add_node(&Vec2::new(150.0, 150.0), ());
        graph.add_edge(0, 1, ());
        graph.add_edge(2, 3, ());

        for edge in graph.edges.iter() {
            println!("{:?}", edge);
        }
    }
}

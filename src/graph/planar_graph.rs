use std::collections::HashMap;

use crate::geometry::utils;

use crate::graph::planar_graph::PlanarGraphError::NoSuchVertex;
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

    pub fn get_node_indices(&self) -> (NodeIndex, NodeIndex) {
        (self.src, self.dst)
    }

    pub fn get_data(&self) -> &T {
        &self.data
    }

    pub fn get_label(&self) -> &String {
        &self.label
    }

    pub fn includes(&self, id: NodeIndex) -> bool {
        self.src == id || self.dst == id
    }
}

pub enum PlanarGraphError {
    NoSuchVertex,
    NoSuchEdge,
}

/// A planar graph data structure with generic node data `T` and edge data `U`.
#[derive(Clone)]
pub struct PlanarGraph<T, U> {
    /// The nodes (or vertices) of the graph
    nodes: Vec<Node<T>>,

    /// The edges of the graph
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

    pub fn has_node_at_index(&self, id: NodeIndex) -> bool {
        id > 0 && id < self.nodes.len()
    }

    pub fn has_edge_at_index(&self, id: EdgeIndex) -> bool {
        id > 0 && id < self.edges.len()
    }

    /// Returns `true` if an edge exists between 2 nodes with indices `a` and `b` and `false` otherwise.
    pub fn has_edge_between_nodes(&self, a: NodeIndex, b: NodeIndex) -> bool {
        for edge in self.edges.iter() {
            // Order does not matter in this graph implementation (it is undirected)
            if (a, b) == (edge.src, edge.dst) || (b, a) == (edge.src, edge.dst) {
                return true;
            }
        }
        false
    }

    /// Returns the number of edges that touch (i.e. include) the specified node, or `None` if the node
    /// does not exist.
    pub fn number_of_edges_incident_to(&self, id: NodeIndex) -> Option<usize> {
        if let Some(node) = self.get_vertex(id) {
            return Some(self.edges.iter().filter(|edge| edge.includes(id)).count());
        }
        None
    }

    /// Maps the graph's nodes into a list of 2D positions.
    fn gather_node_positions(&self) -> Vec<Vec2> {
        self.nodes
            .iter()
            .map(|node| node.position)
            .collect::<Vec<_>>()
    }

    /// Attempts to add a new node to the graph at `position`. This function returns a tuple with 2 entries:
    ///
    /// 1. The index of the newly added node (or the index of an existing node if the specified node was
    ///    close to an existing node)
    /// 2. A list containing the indices of all of the edges that changed as a result of adding the new node
    pub fn add_node(&mut self, position: &Vec2, data: T) -> (NodeIndex, Vec<EdgeIndex>) {
        // Check if the vertex is the same as an existing vertex (within epsilon)
        let (index, distance) = utils::find_closest_to(position, &self.gather_node_positions());

        if distance < nalgebra_glm::epsilon() {
            // If the found distance is less than the specified threshold, the specified
            // node is considered a "duplicate," so we return the index of the existing
            // node and an empty array (since by definition, no edges were modified)
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
        if !self.has_edge_between_nodes(src, dst) {
            // Make sure `src` and `dst` vertices actually exist
            if let (Some(_), Some(_)) = (self.get_vertex(src), self.get_vertex(dst)) {
                // First, push back the new edge
                self.edges.push(Edge::new(src, dst, data, DEFAULT_LABEL));

                // TODO: if the edge we just added gets split, do we copy its data over to the child edges?

                // Perform line-segment/line-segment intersection tests
                let (changed_nodes, mut changed_edges) =
                    self.split_edges_along_edge(self.edge_count() - 1);

                // If the array of changed edges is empty, this means that no edge splitting was necessary,
                // but we still want to make sure we return at least one edge index (in this case, the edge
                // was simply added to the end of the graph's edge array, so just return that index)
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
        if let Some(_) = self.nodes.get(id) {
            println!("Splitting edges at node with index {}", id);
            let mut changed_edges = vec![];
            let mut additional_edges = vec![];

            for (index, edge) in self.edges.iter_mut().enumerate() {
                // Does the target node lie along this edge somewhere?
                let on_edge = utils::is_on_line_segment(
                    &self.nodes[edge.src].position,
                    &self.nodes[edge.dst].position,
                    &self.nodes[id].position,
                    false,
                );

                if on_edge {
                    println!("Node splits edge with index {}", index);

                    // This edge isn't added to the graph's list of edges right away for 2 reasons:
                    //
                    // 1. Rust's borrow checker
                    // 2. We don't need to (re)test if the target node is on this edge, since we
                    //    just created it
                    //
                    // Note that we copy the original edge's data to the 2 new edges (here and below),
                    // but we do NOT copy its label, which would be invalid at this point
                    additional_edges.push(Edge::new(edge.dst, id, *edge.get_data(), DEFAULT_LABEL));

                    // Put one of the new edges in the slot that was previously occupied by the old,
                    // un-split edge and the other at the end of the array - the prior is done so that
                    // we don't have to worry about the rest of the edges being "shuffled" as a result
                    // of a standard array "remove" operation...essentially, we want to add/remove edges
                    // in-place
                    *edge = Edge::new(edge.src, id, *edge.get_data(), DEFAULT_LABEL);
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

            return changed_edges;
        }
        vec![]
    }

    /// Splits any edges that intersect with the specified edge. This function returns a tuple of 2 arrays:
    /// 1. An array containing the indices of any newly created nodes
    /// 2. An array containing the indices of any newly created (or modified / invalidated) edges
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
                    &self.nodes[edge.get_node_indices().0].position,
                    &self.nodes[edge.get_node_indices().1].position,
                    &self.nodes[self.edges[id].get_node_indices().0].position,
                    &self.nodes[self.edges[id].get_node_indices().1].position,
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
            let (index, distance) =
                utils::find_closest_to(intersection, &self.gather_node_positions());

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

    /// Removes an edge from the planar graph, potentially removing any stray nodes as well.
    pub fn remove_edge(&mut self, id: EdgeIndex) -> (Vec<NodeIndex>, Vec<EdgeIndex>) {
        // Does this edge actually exist?
        if let Some(edge) = self.edges.get(id) {
            let (index_a, index_b) = edge.get_node_indices();

            // If either of this edge's endpoints would be stray nodes upon this edge's deletion,
            // we can simply re-use the node deletion procedure below, which will remove this edge
            // as a side effect
            if self.number_of_edges_incident_to(index_a).unwrap() == 1 {
                return self.remove_node(index_a);
            } else if self.number_of_edges_incident_to(index_b).unwrap() == 1 {
                return self.remove_node(index_b);
            }

            // Otherwise, both of the endpoints of this edge are also part of some other edges, so
            // we can simply delete the edge itself
            let mut changed_edges = vec![id];

            // If there is at least one other edge besides this one, replace the edge to-be-deleted
            // with the last edge and return both indices - doing this prevents the entire array from
            // being shuffled
            //
            // Note that the case where there is *only* 1 edge in the graph should be handled by the
            // previous logic (i.e. the calls to `remove_node`)
            // TODO: this could be simplified - do we even need this if-statement here? aren't we guaranteed
            //  to have more than 1 edge at this point?
            if self.edge_count() > 1 {
                // Make sure to do this BEFORE the call to `pop()`, which changes the array length
                changed_edges.push(self.edge_count() - 1);

                let last = self.edges.pop().expect("This should never happen");
                self.edges[id] = last;
            }

            // No nodes should ever be affected by this procedure if we've reached this point
            return (vec![], changed_edges);
        }

        (vec![], vec![])
    }

    /// Removes a node and any incident edges from the planar graph, potentially removing any extra stray
    /// nodes as well.
    pub fn remove_node(&mut self, id: NodeIndex) -> (Vec<NodeIndex>, Vec<EdgeIndex>) {
        // Does this node actually exist?
        if let Some(_) = self.nodes.get(id) {
            // This list will contain the indices of all of the nodes that need to be deleted
            let mut marked_nodes = vec![id];

            // Remove any edges that point to the deleted node
            let (stray_nodes, marked_edges) = self.remove_edges_incident_to_node(id);
            marked_nodes.extend(&stray_nodes);

            // A dictionary from "old" to "new" node indices - the node in the 5th position of the
            // array, for example, will shift downwards some amount if we delete nodes before it, and
            // this data structure captures those relationships for all nodes that are still present
            // after the deletion operation
            let mut remapped_nodes = HashMap::new();

            for (index, _) in self.nodes.iter().enumerate() {
                // If this node was not marked for deletion, its new index in the array will be its
                // current position minus the number of to-be-deleted nodes that come *before* it in
                // the array
                if !marked_nodes.contains(&index) {
                    let back = marked_nodes.iter().filter(|&entry| *entry < index).count();

                    remapped_nodes.insert(index, index - back);
                }
            }

            // Keep track of the index of the "leftmost" node / edge to-be-deleted
            let min_node_index = marked_nodes.iter().cloned().min().unwrap();
            let min_edge_index = marked_edges.iter().cloned().min().unwrap();

            // All nodes / edges that are to the "right" of the indices calculated above will change
            let changed_nodes = (min_node_index..self.node_count()).collect::<Vec<_>>();
            let changed_edges = (min_edge_index..self.edge_count()).collect::<Vec<_>>();

            // Perform the actual deletion operation: the code below deletes all of the marked nodes (or edges)
            // in-place (i.e. simultaneously)
            self.nodes = self
                .nodes
                .iter()
                .cloned()
                .enumerate()
                .filter(|(index, _)| !marked_nodes.contains(index)) // If this node index isn't marked for deletion, keep it
                .map(|index_and_node| index_and_node.1) // Get the second entry of each tuple (just the node itself, not its index)
                .collect::<Vec<_>>();

            self.edges = self
                .edges
                .iter()
                .cloned()
                .enumerate()
                .filter(|(index, _)| !marked_edges.contains(index))
                .map(|index_and_edge| index_and_edge.1)
                .collect::<Vec<_>>();

            // Update any edges that pointed to nodes that were shuffled / moved as a result
            // of the deletion procedure
            for edge in self.edges.iter_mut() {
                // If the dictionary of "old" to "new" node indices contains either of this edge's
                // endpoints, update those indices to reflect the new graph structure
                if remapped_nodes.contains_key(&edge.src) {
                    edge.src = remapped_nodes[&edge.src];
                }
                if remapped_nodes.contains_key(&edge.dst) {
                    edge.dst = remapped_nodes[&edge.dst];
                }
            }

            return (changed_nodes, changed_edges);
        }

        (vec![], vec![])
    }

    /// Finds the indices of all of the edges (and stray nodes) that should be marked for deletion after
    /// removing the specified node - note that this function doesn't actually perform the deletion of
    /// any of these objects. This function returns a tuple of 2 arrays:
    /// 1. An array containing the indices of any stray nodes that should be marked for deletion
    /// 2. An array containing the indices of any edges that should be marked for deletion
    fn remove_edges_incident_to_node(&mut self, id: NodeIndex) -> (Vec<NodeIndex>, Vec<EdgeIndex>) {
        let mut marked_nodes = vec![];
        let mut marked_edges = vec![];

        for (index, edge) in self.edges.iter().enumerate() {
            if edge.includes(id) {
                // This edge should be marked for deletion, as it contains the node that we want to delete
                marked_edges.push(index);

                // Deleting an edge may result in a "floating" stray node, which needs to be deleted as well -
                // this is the edge's other node (i.e. not the one that is already marked for deletion)
                let neighbor = if edge.src == id { edge.dst } else { edge.src };
                // TODO: the JS code below was a lot cleaner
                //let neighbor = edge[1 - edge.indexOf(targetIndex)];

                // Is there another edge (other than this one) that includes the specified node? If not, it
                // is a "stray" node
                if let Some(1) = self.number_of_edges_incident_to(neighbor) {
                    marked_nodes.push(neighbor);
                }
            }
        }

        (marked_nodes, marked_edges)
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

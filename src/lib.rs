mod geometry;
mod graph;
mod math;
mod origami;

use graph::graph::PlanarGraph;

use nalgebra_glm::Vec2;
use wasm_bindgen::prelude::*;
use wasm_bindgen::JsCast;

#[wasm_bindgen]
pub fn run() {
    let window = web_sys::window().unwrap();
    let document = window.document().unwrap();
    let body = document.body().unwrap();

    start(&document, &body);
}

fn start(document: &web_sys::Document, body: &web_sys::HtmlElement) -> Result<(), JsValue> {
    // Create some text elements
    let title = document.create_element("h1")?;
    title.set_inner_html("かみ");
    let clicked = document.create_element("p")?;
    clicked.set_inner_html("No node selected");
    clicked.set_id("clicked");

    // Append the text elements to the body
    body.append_child(&title)?;
    body.append_child(&clicked)?;

    let mut graph: PlanarGraph<(), ()> = PlanarGraph::new();
    graph.add_node(&Vec2::new(50.0, 100.0), ());
    graph.add_node(&Vec2::new(50.0, 200.0), ());
    graph.add_node(&Vec2::new(50.0, 300.0), ());

    graph.add_edge(0, 1, ());

    // Create SVG elements
    let svg = document
        .create_element_ns(Some("http://www.w3.org/2000/svg"), "svg")?
        .dyn_into::<web_sys::SvgElement>()?;

    svg.style().set_property("border", "solid")?;
    svg.set_attribute("width", "500")?;
    svg.set_attribute("height", "500")?;
    svg.set_attribute("viewBox", "0 0 500 500")?;

    // Draw nodes
    for (i, node) in graph.get_nodes().iter().enumerate() {
        // See: https://rustwasm.github.io/wasm-bindgen/examples/closures.html
        let p = document
            .get_element_by_id("clicked")
            .expect("Should have an element with ID #clicked somewhere on the page");

        let callback = Closure::wrap(Box::new(move || {
            p.set_inner_html(&format!("Node {}", i));
        }) as Box<dyn FnMut()>);

        let element = document.create_element_ns(Some("http://www.w3.org/2000/svg"), "circle")?;
        element.set_attribute("cx", &node.get_position().x.to_string())?;
        element.set_attribute("cy", &node.get_position().y.to_string())?;
        element.set_attribute("r", "20")?;
        element.set_attribute("stroke", "black")?;
        element.set_attribute("fill", "blue")?;
        element
            .dyn_ref::<web_sys::SvgElement>()
            .expect("Could not cast to SvgElement type")
            .set_onclick(Some(callback.as_ref().unchecked_ref()));

        svg.append_child(&element)?;

        // See reference code above
        callback.forget();
    }

    // Draw edges
    for edge in graph.get_edges().iter() {
        let (src, dst) = edge.get_indices();
        let src_vertex = graph.get_vertex(src).unwrap();
        let dst_vertex = graph.get_vertex(dst).unwrap();

        let element = document.create_element_ns(Some("http://www.w3.org/2000/svg"), "line")?;
        element.set_attribute("x1", &src_vertex.get_position().x.to_string())?;
        element.set_attribute("y1", &src_vertex.get_position().y.to_string())?;
        element.set_attribute("x2", &dst_vertex.get_position().x.to_string())?;
        element.set_attribute("y2", &dst_vertex.get_position().y.to_string())?;
        element.set_attribute("stroke", "black")?;

        svg.append_child(&element)?;
    }

    body.append_child(&svg)?;
    Ok(())
}

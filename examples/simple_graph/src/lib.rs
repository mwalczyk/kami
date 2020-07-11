use kami::glm::Vec2;
use kami::graph::planar_graph::PlanarGraph;

use wasm_bindgen::prelude::*;
use wasm_bindgen::JsCast;

#[wasm_bindgen]
pub fn run() -> Result<(), JsValue> {
    // Grab references to the HTML document, etc.
    let window = web_sys::window().unwrap();
    let document = window.document().unwrap();
    let body = document.body().unwrap();

    // Create the title
    let heading = document.create_element("h1")?;
    heading.set_inner_html("かみ");
    body.append_child(&heading)?;

    // Create the main SVG element
    let svg = document
        .create_element_ns(Some("http://www.w3.org/2000/svg"), "svg")?
        .dyn_into::<web_sys::SvgElement>()?;
    svg.style().set_property("border", "solid")?;
    svg.set_attribute("width", "500")?;
    svg.set_attribute("height", "500")?;
    svg.set_attribute("viewBox", "0 0 500 500")?;
    body.append_child(&svg)?;

    // Per-example drawing / additions
    draw(&document, &body, &svg)?;

    Ok(())
}

pub fn draw(
    document: &web_sys::Document,
    body: &web_sys::HtmlElement,
    svg: &web_sys::SvgElement,
) -> Result<(), JsValue> {
    // Construct a simple graph
    let mut graph: PlanarGraph<(), ()> = PlanarGraph::new();
    graph.add_node(&Vec2::new(50.0, 100.0), ());
    graph.add_node(&Vec2::new(50.0, 200.0), ());
    graph.add_node(&Vec2::new(50.0, 300.0), ());

    graph.add_edge(0, 1, ());

    // Append a new text element to the body that will tell us which node we have clicked
    let clicked = document.create_element("p")?;
    clicked.set_inner_html("No node selected");
    clicked.set_id("clicked");
    body.append_child(&clicked)?;

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
        element.set_attribute("fill", "red")?;
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

    Ok(())
}

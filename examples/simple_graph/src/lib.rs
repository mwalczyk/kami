mod shared;

use kami::glm::Vec2;
use kami::graph::planar_graph::PlanarGraph;

use wasm_bindgen::prelude::*;
use wasm_bindgen::JsCast;

#[wasm_bindgen]
pub fn run() -> Result<(), JsValue> {
    shared::set_panic_hook();

    // Grab references to the HTML document, etc.
    let window = web_sys::window().unwrap();
    let document = window.document().unwrap();
    let body = document.body().unwrap();

    // Create the title
    let heading = document.create_element("h1")?;
    heading.set_inner_html("かみ");
    body.append_child(&heading)?;

    // Create the main SVG element
    let w = 500.0f32;
    let h = 500.0f32;
    let svg = shared::svg_main(&document, w, h).unwrap();
    body.append_child(&svg)?;

    // Per-example drawing / additions
    draw(&document, &body, &svg, w, h)?;

    web_sys::console::log_1(&"Finished drawing".into());

    Ok(())
}

pub fn draw(
    document: &web_sys::Document,
    body: &web_sys::HtmlElement,
    svg: &web_sys::SvgElement,
    w: f32,
    h: f32,
) -> Result<(), JsValue> {
    // Construct a simple graph
    let mut graph: PlanarGraph<(), ()> = PlanarGraph::new();
    graph.add_node(&Vec2::new(100.0, 100.0), ());
    graph.add_node(&Vec2::new(200.0, 100.0), ());
    graph.add_node(&Vec2::new(150.0, 50.0), ());
    graph.add_node(&Vec2::new(150.0, 150.0), ());
    graph.add_edge(0, 1, ());
    graph.add_edge(2, 3, ());

    // Create a new text element that tells us some basic info about this graph
    let graph_stats = document.create_element("p")?;
    graph_stats.set_inner_html(&format!(
        "Graph has {} edges and {} nodes",
        graph.edge_count(),
        graph.node_count()
    ));
    body.append_child(&graph_stats)?;

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

        let element = shared::svg_circle(&document, &node.get_position(), 5.0, "green", "black").unwrap();
        element
            .dyn_ref::<web_sys::SvgElement>()
            .expect("Could not cast to SvgElement type")
            .set_onclick(Some(callback.as_ref().unchecked_ref()));

        let title = shared::svg_title(&document, &format!("Node {}", i)).unwrap();
        element.append_child(&title);

        svg.append_child(&element)?;

        // See reference code above
        callback.forget();
    }

    // Draw edges
    for edge in graph.get_edges().iter() {
        let (src, dst) = edge.get_indices();
        let src_vertex = graph.get_vertex(src).unwrap();
        let dst_vertex = graph.get_vertex(dst).unwrap();

        let element = shared::svg_line(
            &document,
            src_vertex.get_position(),
            dst_vertex.get_position(),
            "black",
        )
        .unwrap();

        svg.append_child(&element)?;
    }

    Ok(())
}

use kami::geometry::utils;
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
    let w = 500.0f32;
    let h = 500.0f32;
    let svg = document
        .create_element_ns(Some("http://www.w3.org/2000/svg"), "svg")?
        .dyn_into::<web_sys::SvgElement>()?;
    svg.style().set_property("border", "solid")?;
    svg.set_attribute("width", &w.to_string())?;
    svg.set_attribute("height", &h.to_string())?;
    svg.set_attribute("viewBox", "0 0 500 500")?;
    body.append_child(&svg)?;

    // Per-example drawing / additions
    draw(&document, &body, &svg, w, h)?;

    web_sys::console::log_1(&"Finished drawing".into());

    Ok(())
}

pub fn svg_polygon(
    document: &web_sys::Document,
    points: &Vec<Vec2>,
    fill: &str,
    stroke: &str,
) -> Result<web_sys::Element, JsValue> {
    // Convert vec of points into a single string like: "200,300 100,200 400,50"
    let mut points_attr = points
        .iter()
        .map(|point| format!("{},{} ", point.x, point.y))
        .collect::<String>();

    // Remove the last space (probably doesn't matter, but we do it anyways)
    points_attr.pop();

    let element = document.create_element_ns(Some("http://www.w3.org/2000/svg"), "polygon")?;
    element.set_attribute("points", &points_attr)?;
    element.set_attribute("fill", fill)?;
    element.set_attribute("stroke", stroke)?;

    Ok(element)
}

pub fn svg_circle(
    document: &web_sys::Document,
    position: &Vec2,
    radius: f32,
    fill: &str,
    stroke: &str,
) -> Result<web_sys::Element, JsValue> {
    let element = document.create_element_ns(Some("http://www.w3.org/2000/svg"), "circle")?;
    element.set_attribute("cx", &position.x.to_string())?;
    element.set_attribute("cy", &position.y.to_string())?;
    element.set_attribute("r", &radius.to_string())?;
    element.set_attribute("fill", fill)?;
    element.set_attribute("stroke", stroke)?;

    Ok(element)
}

pub fn svg_line(
    document: &web_sys::Document,
    a: &Vec2,
    b: &Vec2,
    stroke: &str,
) -> Result<web_sys::Element, JsValue> {
    let element = document.create_element_ns(Some("http://www.w3.org/2000/svg"), "line")?;
    element.set_attribute("x1", &a.x.to_string())?;
    element.set_attribute("y1", &a.y.to_string())?;
    element.set_attribute("x2", &b.x.to_string())?;
    element.set_attribute("y2", &b.y.to_string())?;
    element.set_attribute("stroke", stroke)?;

    Ok(element)
}

pub fn draw(
    document: &web_sys::Document,
    body: &web_sys::HtmlElement,
    svg: &web_sys::SvgElement,
    w: f32,
    h: f32,
) -> Result<(), JsValue> {
    // A simple triangle, 3 points
    let points = vec![
        Vec2::new(w * 0.5 + 200.0, h * 0.5 - 50.0),
        Vec2::new(w * 0.5 - 200.0, h * 0.5 - 50.0),
        Vec2::new(w * 0.5, h * 0.5 + 150.0),
    ];

    // If debugging...
    //web_sys::console::log_1(&points_attr.into());

    let element = svg_polygon(&document, &points, "white", "black").unwrap();
    svg.append_child(&element)?;

    // Calculate the triangle incenter
    if let Some((incenter, inradius)) =
        utils::calculate_triangle_incenter(&points[0], &points[1], &points[2])
    {
        // Draw the incenter
        let element = svg_circle(&document, &incenter, inradius, "white", "blue").unwrap();
        svg.append_child(&element)?;

        // Draw the angle bisectors, which (by definition) meet at the triangle incenter
        let line_a = svg_line(&document, &points[0], &incenter, "red").unwrap();
        let line_b = svg_line(&document, &points[1], &incenter, "red").unwrap();
        let line_c = svg_line(&document, &points[2], &incenter, "red").unwrap();
        svg.append_child(&line_a)?;
        svg.append_child(&line_b)?;
        svg.append_child(&line_c)?;
    }

    Ok(())
}

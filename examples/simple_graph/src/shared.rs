use kami::glm::Vec2;

use wasm_bindgen::prelude::*;
use wasm_bindgen::JsCast;

pub fn set_panic_hook() {
    // https://github.com/rustwasm/console_error_panic_hook#readme
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

pub fn svg_main(document: &web_sys::Document, w: f32, h: f32) -> Result<web_sys::SvgElement, JsValue> {
    let svg = document
        .create_element_ns(Some("http://www.w3.org/2000/svg"), "svg")?
        .dyn_into::<web_sys::SvgElement>()?;
    svg.style().set_property("border", "solid")?;
    svg.set_attribute("width", &w.to_string())?;
    svg.set_attribute("height", &h.to_string())?;
    svg.set_attribute("viewBox", &format!("0 0 {} {}", w as i32, h as i32))?;

    Ok(svg)
}

pub fn svg_title(document: &web_sys::Document, text: &str) -> Result<web_sys::Element, JsValue> {
    let element = document.create_element_ns(Some("http://www.w3.org/2000/svg"), "title")?;
    element.set_inner_html(text);

    Ok(element)
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
mod shared;

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

    Ok(())
}

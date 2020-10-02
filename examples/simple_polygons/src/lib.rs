mod shared;

use kami::geometry::polygon::Polygon;
use kami::geometry::utils;
use kami::glm::Vec2;
use kami::graph::planar_graph::PlanarGraph;

use rand::Rng;
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
    // PRNG
    let mut rng = rand::thread_rng();

    // Create some random triangles, 3 points each
    let mut triangles = vec![];
    let num_triangles = 5;
    for _ in 0..num_triangles {
        let a = Vec2::new(rng.gen_range(0.0, w * 0.5), rng.gen_range(h * 0.5, h));
        let b = Vec2::new(rng.gen_range(0.0, w * 0.5), rng.gen_range(h * 0.5, h));
        let c = Vec2::new(rng.gen_range(0.0, w * 0.5), rng.gen_range(h * 0.5, h));
        triangles.push((a, b, c));
    }

    // Draw each triangle and its incenter
    for (a, b, c) in triangles.iter() {
        let element = shared::svg_polygon(&document, &vec![*a, *b, *c], "none", "black").unwrap();
        svg.append_child(&element)?;

        // Calculate the triangle incenter
        if let Some((incenter, inradius)) = utils::calculate_triangle_incenter(a, b, c) {
            // Draw the incenter
            let element =
                shared::svg_circle(&document, &incenter, inradius, "none", "blue").unwrap();
            svg.append_child(&element)?;

            // Draw the angle bisectors, which (by definition) meet at the triangle incenter
            let line_a = shared::svg_line(&document, a, &incenter, "red").unwrap();
            let line_b = shared::svg_line(&document, b, &incenter, "red").unwrap();
            let line_c = shared::svg_line(&document, c, &incenter, "red").unwrap();
            svg.append_child(&line_a)?;
            svg.append_child(&line_b)?;
            svg.append_child(&line_c)?;
        }
    }
    // If debugging...
    //web_sys::console::log_1(&points_attr.into());

    // Create some random line segments
    let mut segments = vec![];
    let num_segments = 100;
    for _ in 0..num_segments {
        let a = Vec2::new(rng.gen_range(0.0, w), rng.gen_range(0.0, h * 0.5));
        let b = Vec2::new(rng.gen_range(0.0, w), rng.gen_range(0.0, h * 0.5));
        segments.push((a, b));
    }

    // Find intersections between segments
    let mut intersections = 0;
    for i in 0..segments.len() {
        // Create the i-th line segment
        let element = shared::svg_line(&document, &segments[i].0, &segments[i].1, "green").unwrap();
        svg.append_child(&element)?;

        // Loop over all other line segments and check for intersections
        for j in (i + 1)..segments.len() {
            if let Some(intersection) = utils::calculate_line_segment_intersection(
                &segments[i].0,
                &segments[i].1,
                &segments[j].0,
                &segments[j].1,
            ) {
                let element = shared::svg_circle(&document, &intersection, 2.0, "red", "").unwrap();
                svg.append_child(&element)?;

                intersections += 1;//
            }
        }
    }

    let count = document.create_element("p")?;
    count.set_inner_html(&format!("Found {} intersections", intersections));
    body.append_child(&count)?;

    // Test triangulation
    let mut polygon = Polygon::new(&vec![
        Vec2::new(0.0, 0.0),
        Vec2::new(10.0, 7.0),
        Vec2::new(12.0, 3.0),
        Vec2::new(20.0, 8.0),
        Vec2::new(13.0, 17.0),
        Vec2::new(10.0, 12.0),
        Vec2::new(12.0, 14.0),
        Vec2::new(14.0, 9.0),
        Vec2::new(8.0, 10.0),
        Vec2::new(6.0, 14.0),
        Vec2::new(10.0, 15.0),
        Vec2::new(7.0, 18.0),
        Vec2::new(0.0, 16.0),
        Vec2::new(1.0, 13.0),
        Vec2::new(3.0, 15.0),
        Vec2::new(5.0, 8.0),
        Vec2::new(-2.0, 9.0),
        Vec2::new(5.0, 5.0),
    ])
    .unwrap();
    //let mut polygon = Polygon::regular(8, 10.0, &Vec2::new(0.0, 0.0));
    polygon.uniform_scale(10.0);
    polygon.translate(&Vec2::new(250.0, 250.0));

    let element = shared::svg_polygon(&document, polygon.get_vertices(), "none", "black").unwrap();
    svg.append_child(&element)?;

    let triangles = polygon.triangulate();

    for (index, triangle) in triangles.iter().enumerate() {
        let element = shared::svg_polygon(
            &document,
            triangle.get_vertices(),
            &format!("hsl({}, 100%, 50%)", (index as f32 / triangles.len() as f32) * 360.0),
            "none",
        )
        .unwrap();
        svg.append_child(&element)?;
    }

    Ok(())
}

use itertools::Itertools;

use rukki::graph_algos::superbubble;
use rukki::*;

#[test]
fn multi_link_bubble() {
    let s = "
S a * LN:i:100
S b * LN:i:100
L a + b + 50M
L a + b + 75M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let bubble = superbubble::find_superbubble(
        &g,
        Vertex::forward(0),
        &superbubble::SbSearchParams::unrestricted(),
    )
    .unwrap();
    assert_eq!(bubble.length_range(&g), (125, 150));
    assert!(g.name(bubble.end_vertex().node_id) == "b");
}

#[test]
#[should_panic]
fn extra_link_start() {
    let s = "
S a * LN:i:100
S c * LN:i:100
S b * LN:i:100
L a + b + 50M
L a + b + 75M
L a + c + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let _bubble = superbubble::find_superbubble(
        &g,
        Vertex::forward(0),
        &superbubble::SbSearchParams::unrestricted(),
    )
    .unwrap();
}

#[test]
#[should_panic]
fn extra_link_end() {
    let s = "
S a * LN:i:100
S b * LN:i:100
S c * LN:i:100
L a + b + 50M
L a + b + 50M
L c + b + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let _bubble = superbubble::find_superbubble(
        &g,
        Vertex::forward(0),
        &superbubble::SbSearchParams::unrestricted(),
    )
    .unwrap();
}

#[test]
fn simple_bubble() {
    let s = "
S a * LN:i:100
S b * LN:i:100
S c * LN:i:100
S d * LN:i:100
L a + b + 50M
L a + c + 50M
L b + d + 50M
L c + d + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let bubble = superbubble::find_superbubble(
        &g,
        Vertex::forward(0),
        &superbubble::SbSearchParams::unrestricted(),
    )
    .unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "d");
    assert_eq!(bubble.vertices().count(), 4);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect_vec();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a+", "b+", "c+", "d+"]);
    assert_eq!(bubble.length_range(&g), (200, 200));
}

//TODO support this case
#[test]
#[should_panic]
fn simple_bubble_loop() {
    let s = "
S a * LN:i:100
S b * LN:i:100
S c * LN:i:100
L a + b + 50M
L a + c + 50M
L b + a + 50M
L c + a + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let bubble = superbubble::find_superbubble(
        &g,
        Vertex::forward(0),
        &superbubble::SbSearchParams::unrestricted(),
    )
    .unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "a");
    assert_eq!(bubble.vertices().count(), 4);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect_vec();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a+", "b+", "c+"]);
    assert_eq!(bubble.length_range(&g), (100, 100));
}

#[test]
fn triple_bubble() {
    let s = "
S a * LN:i:100
S b * LN:i:100
S c * LN:i:100
S d * LN:i:100
L a + b + 50M
L a + c + 50M
L b + d + 50M
L c + d + 50M
L a + d + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let bubble = superbubble::find_superbubble(
        &g,
        Vertex::forward(0),
        &superbubble::SbSearchParams::unrestricted(),
    )
    .unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "d");
    assert_eq!(bubble.vertices().count(), 4);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect_vec();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a+", "b+", "c+", "d+"]);
    assert_eq!(bubble.length_range(&g), (150, 200));
}

#[test]
fn super_bubble_1() {
    let s = "
S a * LN:i:100
S b * LN:i:100
S c * LN:i:100
S d * LN:i:100
L a + b + 50M
L a + c + 50M
L b + c + 50M
L b + d + 50M
L c + d + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let bubble = superbubble::find_superbubble(
        &g,
        Vertex::forward(0),
        &superbubble::SbSearchParams::unrestricted(),
    )
    .unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "d");
    assert_eq!(bubble.vertices().count(), 4);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect_vec();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a+", "b+", "c+", "d+"]);
    assert_eq!(bubble.length_range(&g), (200, 250));
}

#[test]
fn super_bubble_1_reverse() {
    let s = "
S a * LN:i:100
S b * LN:i:100
S c * LN:i:100
S d * LN:i:100
L a + b + 50M
L a + c + 50M
L b + c + 50M
L b + d + 50M
L c + d + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let bubble = superbubble::find_superbubble(
        &g,
        Vertex::reverse(3),
        &superbubble::SbSearchParams::unrestricted(),
    )
    .unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "a");
    assert_eq!(bubble.vertices().count(), 4);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect_vec();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a-", "b-", "c-", "d-"]);
    assert_eq!(bubble.length_range(&g), (200, 250));
}

#[test]
fn super_bubble_2() {
    let s = "
S a * LN:i:100
S b * LN:i:100
S c * LN:i:100
S d * LN:i:100
S e * LN:i:100
S f * LN:i:100
L a + b + 75M
L a + c + 50M
L b + d + 75M
L b + e + 50M
L c + d + 50M
L c + e + 50M
L d + f + 50M
L e + f + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let bubble = superbubble::find_superbubble(
        &g,
        Vertex::forward(0),
        &superbubble::SbSearchParams::unrestricted(),
    )
    .unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "f");
    assert_eq!(bubble.vertices().count(), 6);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect_vec();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a+", "b+", "c+", "d+", "e+", "f+"]);
    assert_eq!(bubble.length_range(&g), (200, 250));
}

#[test]
fn simple_max_chain() {
    let s = "
S a * LN:i:100
S b * LN:i:100
S c * LN:i:100
S d * LN:i:100
S e * LN:i:100
S f * LN:i:100
S g * LN:i:100
L a + b + 50M
L a + c + 50M
L b + d + 50M
L c + d + 50M
L d + e + 50M
L d + f + 50M
L e + g + 50M
L f + g + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    let chain = superbubble::find_maximal_chain(
        &g,
        Vertex::forward(g.name2id("d")),
        &superbubble::SbSearchParams::unrestricted(),
    );
    assert_eq!(chain.len(), 2);
    assert_eq!(chain[0].start_vertex(), Vertex::forward(g.name2id("a")));
    assert_eq!(chain[0].end_vertex(), Vertex::forward(g.name2id("d")));
    assert_eq!(chain[1].start_vertex(), Vertex::forward(g.name2id("d")));
    assert_eq!(chain[1].end_vertex(), Vertex::forward(g.name2id("g")));
    assert_eq!(superbubble::length_range(&chain, &g), (300, 300));
}

#[test]
fn simple_chain_loop() {
    let s = "
S a * LN:i:100
S b * LN:i:100
S c * LN:i:100
S d * LN:i:100
S e * LN:i:100
S f * LN:i:100
L a + b + 50M
L a + c + 50M
L b + d + 50M
L c + d + 50M
L d + e + 50M
L d + f + 50M
L e + a + 50M
L f + a + 50M
";
    let g = Graph::read(&s.replace(' ', "\t"));
    //testing search ahead
    let chain = superbubble::find_chain_ahead(
        &g,
        Vertex::forward(0),
        &superbubble::SbSearchParams::unrestricted(),
    );
    assert_eq!(chain.len(), 2);
    assert!(g.name(chain[0].end_vertex().node_id) == "d");
    assert!(g.name(chain[1].end_vertex().node_id) == "a");
    assert_eq!(superbubble::length_range(&chain, &g), (200, 200));

    //testing maximal chain search
    let chain = superbubble::find_maximal_chain(
        &g,
        Vertex::forward(g.name2id("d")),
        &superbubble::SbSearchParams::unrestricted(),
    );
    assert_eq!(chain.len(), 2);
    assert_eq!(chain[0].start_vertex(), Vertex::forward(g.name2id("d")));
    assert_eq!(chain[0].end_vertex(), Vertex::forward(g.name2id("a")));
    assert_eq!(chain[1].start_vertex(), Vertex::forward(g.name2id("a")));
    assert_eq!(chain[1].end_vertex(), Vertex::forward(g.name2id("d")));
    assert_eq!(superbubble::length_range(&chain, &g), (200, 200));
}

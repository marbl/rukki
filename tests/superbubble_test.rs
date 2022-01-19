use graph_analysis::*;
use graph_analysis::graph_algos::superbubble;

#[test]
fn multi_link_bubble() {
    let s = "
S a * LN:i:100
S b * LN:i:100
L a + b + 50M
L a + b + 75M
";
    let g = Graph::read(&s.replace(" ", "\t"));
    let bubble = superbubble::find_superbubble(&g,
        Vertex::forward(0), &superbubble::SbSearchParams::unrestricted()).unwrap();
    assert_eq!(bubble.length_range(), (25, 50));
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
    let g = Graph::read(&s.replace(" ", "\t"));
    let _bubble = superbubble::find_superbubble(&g,
        Vertex::forward(0), &superbubble::SbSearchParams::unrestricted()).unwrap();
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
    let g = Graph::read(&s.replace(" ", "\t"));
    let _bubble = superbubble::find_superbubble(&g,
        Vertex::forward(0), &superbubble::SbSearchParams::unrestricted()).unwrap();
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
    let g = Graph::read(&s.replace(" ", "\t"));
    let bubble = superbubble::find_superbubble(&g,
        Vertex::forward(0), &superbubble::SbSearchParams::unrestricted()).unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "d");
    assert_eq!(bubble.vertices().count(), 4);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect::<Vec<String>>();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a+","b+","c+","d+"]);
    assert_eq!(bubble.length_range(), (100, 100));
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
    let g = Graph::read(&s.replace(" ", "\t"));
    let bubble = superbubble::find_superbubble(&g,
        Vertex::forward(0), &superbubble::SbSearchParams::unrestricted()).unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "d");
    assert_eq!(bubble.vertices().count(), 4);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect::<Vec<String>>();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a+","b+","c+","d+"]);
    assert_eq!(bubble.length_range(), (50, 100));
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
    let g = Graph::read(&s.replace(" ", "\t"));
    let bubble = superbubble::find_superbubble(&g,
        Vertex::forward(0), &superbubble::SbSearchParams::unrestricted()).unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "d");
    assert_eq!(bubble.vertices().count(), 4);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect::<Vec<String>>();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a+","b+","c+","d+"]);
    assert_eq!(bubble.length_range(), (100, 150));
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
    let g = Graph::read(&s.replace(" ", "\t"));
    let bubble = superbubble::find_superbubble(&g,
        Vertex::reverse(3), &superbubble::SbSearchParams::unrestricted()).unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "a");
    assert_eq!(bubble.vertices().count(), 4);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect::<Vec<String>>();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a-","b-","c-","d-"]);
    assert_eq!(bubble.length_range(), (100, 150));
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
    let g = Graph::read(&s.replace(" ", "\t"));
    let bubble = superbubble::find_superbubble(&g,
        Vertex::forward(0), &superbubble::SbSearchParams::unrestricted()).unwrap();
    assert!(g.name(bubble.end_vertex().node_id) == "f");
    assert_eq!(bubble.vertices().count(), 6);
    let mut bubble_vertices = bubble.vertices().map(|&v| g.v_str(v)).collect::<Vec<String>>();
    bubble_vertices.sort();
    assert_eq!(bubble_vertices, vec!["a+","b+","c+","d+","e+","f+"]);
    assert_eq!(bubble.length_range(), (100, 150));
}

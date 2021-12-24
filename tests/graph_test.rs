use graph_analysis::*;

#[test]
fn one_node() {
    let s = "S a * LN:i:100";
    let g = Graph::read(&s.replace(" ", "\t"));
    assert_eq!(1, g.node_cnt());
    assert_eq!(0, g.link_cnt());
    let n = g.all_nodes().next().unwrap();
    assert_eq!("a", n.name);
    assert_eq!(100, n.length);
    assert_eq!(None, g.all_links().next());
}

#[test]
fn loop1() {
    let s = "
S a * LN:i:100
L a + a + 10M
";
    let g = Graph::read(&s.replace(" ", "\t"));
    assert_eq!(1, g.node_cnt());
    assert_eq!(1, g.link_cnt());
    let l = g.all_links().next().unwrap();
    assert_eq!(10, l.overlap);
    assert_eq!(Direction::FORWARD, l.start.direction);
    assert_eq!(Direction::FORWARD, l.end.direction);
}

#[test]
#[should_panic]
fn nontrivial_cigar() {
    let s = "
S a * LN:i:100
L a + a + 1D10M1I
";
    Graph::read(&s.replace(" ", "\t"));
}


#[test]
fn loop2() {
    let s = "
S a * LN:i:100
L a - a - 10M
";
    let g = Graph::read(&s.replace(" ", "\t"));
    assert_eq!(1, g.node_cnt());
    assert_eq!(1, g.link_cnt());
    let l = g.all_links().next().unwrap();
    assert_eq!("a+->a+", g.l_str(l));
    assert_eq!("a", g.node(l.start.node_id).name);
    assert_eq!("a", g.node(l.end.node_id).name);
    assert_eq!(Direction::FORWARD, l.start.direction);
    assert_eq!(Direction::FORWARD, l.end.direction);
}

#[test]
fn self_conj1() {
    let s = "
S a * LN:i:100
L a + a - 10M
";
    let g = Graph::read(&s.replace(" ", "\t"));
    assert_eq!(1, g.node_cnt());
    assert_eq!(1, g.link_cnt());
    let l = g.all_links().next().unwrap();
    assert_eq!("a+->a-", g.l_str(l));
    assert_eq!(Direction::FORWARD, l.start.direction);
    assert_eq!(Direction::REVERSE, l.end.direction);
}

#[test]
fn self_conj2() {
    let s = "
S a * LN:i:100
L a - a + 10M
";
    let g = Graph::read(&s.replace(" ", "\t"));
    assert_eq!(1, g.node_cnt());
    assert_eq!(1, g.link_cnt());
    let l = g.all_links().next().unwrap();
    assert_eq!("a-->a+", g.l_str(l));
    assert_eq!(Direction::REVERSE, l.start.direction);
    assert_eq!(Direction::FORWARD, l.end.direction);
}

#[test]
fn two_nodes() {
    let s = "
S a * LN:i:100
S b * LN:i:200
";
    let g = Graph::read(&s.replace(" ", "\t"));
    assert_eq!(2, g.node_cnt());
    assert_eq!(0, g.link_cnt());
}

#[test]
fn one_link() {
    let s = "
S a * LN:i:100
S b * LN:i:200
L a + b + 10M
";
    let g = Graph::read(&s.replace(" ", "\t"));
    assert_eq!(2, g.node_cnt());
    assert_eq!(1, g.link_cnt());
}

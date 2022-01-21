use std::str;
use std::collections::HashMap;

//TODO which ones are redundant?
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub enum Direction {
    FORWARD,
    REVERSE,
}

impl Direction {
    fn flip(d: Direction) -> Direction {
        match d {
            Self::FORWARD => Self::REVERSE,
            Self::REVERSE => Self::FORWARD,
        }
    }

    fn parse_char(c: char) -> Direction {
        match c {
            '+' => Self::FORWARD,
            '-' => Self::REVERSE,
            _ => panic!("Unknown direction {}", c),
        }
    }

    fn parse(s: &str) -> Direction {
        assert!(s.len() == 1, "Unknown direction {}", s);
        Self::parse_char(s.chars().next().unwrap())
    }

    fn str(d: Direction) -> &'static str {
        match d {
            Self::FORWARD => "+",
            Self::REVERSE => "-",
        }
    }

    fn gaf_str(d: Direction) -> &'static str {
        match d {
            Self::FORWARD => ">",
            Self::REVERSE => "<",
        }
    }
}

#[derive(Clone)]
pub struct Node {
    //node size
    pub name: String,
    pub length: usize,
    pub coverage: f64,
}

//TODO which ones are redundant?
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct Vertex {
    //node id
    pub node_id: usize,
    //direction
    pub direction: Direction,
}

impl Vertex {
    pub fn forward(node_id: usize) -> Vertex {
        Vertex {node_id, direction: Direction::FORWARD}
    }

    pub fn reverse(node_id: usize) -> Vertex {
        Vertex {node_id, direction: Direction::REVERSE}
    }

    pub fn rc(&self) -> Vertex {
        Vertex {
            node_id: self.node_id,
            direction: Direction::flip(self.direction),
        }
    }
}

//FIXME support link coverage!
//TODO separate 'links' and 'edges'
//links will have overlap size, CIGAR, etc
//edges will represent a Vertex pair
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Link {
    pub start: Vertex,
    pub end: Vertex,
    pub overlap: usize,
}

impl Link {
    pub fn rc(&self) -> Link {
        Link {
            start: self.end.rc(),
            end: self.start.rc(),
            overlap: self.overlap,
        }
    }

    fn is_canonical(&self) -> bool {
        self <= &self.rc()
    }

    fn join_same(l1: &Link, l2: &Link) -> bool {
        l1.start == l2.start && l1.end == l2.end
    }

    fn parallel(l1: &Link, l2: &Link) -> bool {
        Self::join_same(l1, l2) || Self::join_same(l1, &l2.rc())
    }
}

pub struct Graph {
    nodes: Vec<Node>,
    //TODO storage is excessive, should only store neighbor
    //incoming & outgoing links for every node
    incoming_links: Vec<Vec<Link>>,
    outgoing_links: Vec<Vec<Link>>,
    //TODO switch to &str and figure out how to work with lifetimes
    name2ids: HashMap<String, usize>,
}

//TODO think about useful iterators and reimplement this one via composition
//FIXME improve when learn how to store iterator as a field :)
struct AllLinkIter<'a> {
    g: &'a Graph,
    curr_node: usize,
    incoming_flag: bool,
    pos: usize,
}

impl<'a> AllLinkIter<'a> {
    fn new(g: &'a Graph) -> AllLinkIter<'a> {
        AllLinkIter {
            g,
            curr_node: 0,
            incoming_flag: true,
            pos: 0,
        }
    }
}

impl<'a> Iterator for AllLinkIter<'a> {
    type Item = Link;

    fn next(&mut self) -> Option<Self::Item> {
        while self.curr_node < self.g.node_cnt() {
            if self.incoming_flag {
                let links = &self.g.incoming_links[self.curr_node];
                assert!(self.pos <= links.len());
                if self.pos < links.len() {
                    let link = links[self.pos];
                    assert!(link.end.node_id == self.curr_node);
                    self.pos += 1;
                    if link.end < link.start {
                        return Some(link);
                    }
                } else {
                    self.incoming_flag = false;
                    self.pos = 0;
                }
            } else {
                let links = &self.g.outgoing_links[self.curr_node];
                assert!(self.pos <= links.len());
                if self.pos < links.len() {
                    let link = links[self.pos];
                    assert!(link.start.node_id == self.curr_node);
                    self.pos += 1;
                    if link.start <= link.end {
                        return Some(link);
                    }
                } else {
                    self.incoming_flag = true;
                    self.pos = 0;
                    self.curr_node += 1;
                }
            }
        }
        return None;
    }
}

struct VertexIter<'a> {
    g: &'a Graph,
    curr_node: usize,
    forward_flag: bool,
}

impl<'a> VertexIter<'a> {
    fn new(g: &'a Graph) -> VertexIter<'a> {
        VertexIter {
            g,
            curr_node: 0,
            forward_flag: true,
        }
    }
}

impl<'a> Iterator for VertexIter<'a> {
    type Item = Vertex;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_node < self.g.node_cnt() {
            if self.forward_flag {
                self.forward_flag = false;
                return Some(Vertex::forward(self.curr_node));
            } else {
                let node_id = self.curr_node;
                self.forward_flag = true;
                self.curr_node += 1;
                return Some(Vertex::reverse(node_id));
            }
        }
        return None;
    }
}

impl Graph {

    pub fn new() -> Graph {
        let g = Graph {
            nodes: Vec::new(),
            incoming_links: Vec::new(),
            outgoing_links: Vec::new(),
            name2ids:  HashMap::new(),
        };

        g
    }

    pub fn node_cnt(&self) -> usize {
        self.nodes.len()
    }

    pub fn node_iter(&self) -> std::slice::Iter<Node> {
        self.nodes.iter()
    }

    pub fn add_node(&mut self, node: Node) -> usize {
        //TODO rewrite without cloning with lifetimes
        let node_id = self.nodes.len();
        self.name2ids.insert(node.name.clone(), node_id);
        self.nodes.push(node);
        self.incoming_links.push(Vec::new());
        self.outgoing_links.push(Vec::new());
        node_id
    }

    pub fn add_link(&mut self, link: Link) {
        //FIXME Currently doesn't check that every link is represented only once
        //TODO Think of some nice 'views' for vectors that will reverse complement everything put
        //there
        match link.start.direction {
            Direction::FORWARD => self.outgoing_links[link.start.node_id].push(link),
            Direction::REVERSE => self.incoming_links[link.start.node_id].push(link.rc()),
        };

        if &link == &link.rc() { return };

        match link.end.direction {
            Direction::FORWARD => self.incoming_links[link.end.node_id].push(link),
            Direction::REVERSE => self.outgoing_links[link.end.node_id].push(link.rc()),
        };
    }

    //FIXME add this check within add_link function
    fn check_links(&self) {
        assert!(self.nodes.len() == self.incoming_links.len());
        assert!(self.nodes.len() == self.outgoing_links.len());
        for (node_id, _) in self.all_nodes().enumerate() {
            let v = Vertex::forward(node_id);
            assert!(self.incoming_links[node_id].iter().filter(|l| l.end != v).count() == 0
                        , "Problem with incoming links for node {}", self.nodes[node_id].name);
            assert!(self.outgoing_links[node_id].iter().filter(|l| l.start != v).count() == 0
                        , "Problem with incoming links for node {}", self.nodes[node_id].name);
        }
    }

    //TODO switch to iterator?
    fn parse_tag<T: str::FromStr>(fields: &[&str], prefix: &str) -> Option<T> {
        fields.iter()
            .filter(|s| s.starts_with(prefix))
            .map(|s| match s[prefix.len()..].parse::<T>() {
                Ok(t) => t,
                Err(_) => panic!("Couldn't parse tag {}", s),
            })
            .next()
    }

    fn parse_overlap(cigar: &str) -> usize {
        assert!(cigar.ends_with('M'), "Invalid overlap {}", cigar);
        let ovl = &cigar[..(cigar.len()-1)];
        ovl.trim().parse().expect(&format!("Invalid overlap {}", cigar))
    }

    //TODO switch to something iterable
    pub fn read(graph_str: &str) -> Graph {
        let mut g = Self::new();

        for line in graph_str.lines() {
            if line.starts_with("S\t") {
                let split: Vec<&str> = line.split('\t').collect();
                //println!("Node line {:?}", split);
                let name = String::from(split[1]);
                let length = if split[2] != "*" {
                                 split[2].trim().len()
                             } else {
                                 Self::parse_tag(&split[3..split.len()], "LN:i:")
                                     .expect("Neither sequence nor LN tag provided")
                             };
                let coverage = match Self::parse_tag::<usize>(&split[3..split.len()], "RC:i:") {
                    None => Self::parse_tag(&split[3..split.len()], "ll:f:")
                                        .unwrap_or(0.),
                    Some(raw_cnt) => raw_cnt as f64 / length as f64,
                };
                g.add_node(Node{name, length, coverage});
            }
        }

        for line in graph_str.lines() {
            if line.starts_with("L\t") {
                let split: Vec<&str> = line.trim().split('\t').collect();
                //println!("Link line {:?}", split);
                let start = Vertex {
                    node_id: g.name2id(split[1]),
                    direction: Direction::parse(split[2]),
                };
                let end = Vertex {
                    node_id: g.name2id(split[3]),
                    direction: Direction::parse(split[4]),
                };
                let overlap = Self::parse_overlap(split[5]);
                g.add_link(Link{start, end, overlap});
            }
        }
        g.check_links();
        g
    }

    pub fn as_gfa(&self) -> String {
        let mut gfa = String::new();

        for n in self.all_nodes() {
            gfa += &format!("S\t{}\t*\tLN:i:{}\tRC:i:{}\tll:f:{:.1}\n",
                n.name, n.length,
                (n.coverage * n.length as f64).round() as u64, n.coverage);
        }

        for l in self.all_links() {
            gfa += &format!("L\t{}\t{}\t{}\t{}\t{}M\n",
                self.node(l.start.node_id).name, Direction::str(l.start.direction),
                self.node(l.end.node_id).name, Direction::str(l.end.direction),
                l.overlap);
        }

        gfa
    }

    //fn get_vertex(&self, name: &str, direction: Direction) -> Vertex {
    //    let node_id = self.name2id(name);
    //    Vertex {node_id, direction}
    //}

    fn rc(links: &Vec<Link>) -> Vec<Link> {
        links.into_iter().map(|x| x.rc()).collect()
    }

    pub fn node(&self, node_id: usize) -> &Node {
        &self.nodes[node_id]
    }

    pub fn node_by_name(&self, name: &str) -> &Node {
        &self.nodes[self.name2id(name)]
    }

    pub fn name(&self, node_id: usize) -> &str {
        &self.node(node_id).name
    }

    //TODO switch to iterators when learn enough Rust :)
    pub fn outgoing_edge_cnt(&self, v: Vertex) -> usize {
        match v.direction {
            Direction::FORWARD => self.outgoing_links[v.node_id].len(),
            Direction::REVERSE => self.incoming_links[v.node_id].len(),
        }
    }

    //TODO switch to iterators when learn enough Rust :)
    pub fn outgoing_edges(&self, v: Vertex) -> Vec<Link> {
        match v.direction {
            Direction::FORWARD => self.outgoing_links[v.node_id].clone(),
            Direction::REVERSE => Self::rc(&self.incoming_links[v.node_id]),
        }
    }

    pub fn incoming_edge_cnt(&self, v: Vertex) -> usize {
        self.outgoing_edge_cnt(v.rc())
    }

    //TODO switch to iterators when learn enough Rust :)
    pub fn incoming_edges(&self, v: Vertex) -> Vec<Link> {
        match v.direction {
            Direction::FORWARD => self.incoming_links[v.node_id].clone(),
            Direction::REVERSE => Self::rc(&self.outgoing_links[v.node_id]),
        }
    }

    pub fn name2id(&self, name: &str) -> usize {
        match self.name2ids.get(name) {
            Some(&id) => id,
            None => panic!("Node {} is not in the graph", name),
        }
    }

    //TODO iterate over references
    pub fn all_links(&self) -> impl Iterator<Item=Link> + '_ {
        AllLinkIter::new(self)
    }

    pub fn all_nodes(&self) -> impl Iterator<Item=&Node> + '_ {
        self.nodes.iter()
    }

    //TODO iterate over references
    pub fn all_vertices(&self) -> impl Iterator<Item=Vertex> + '_ {
        VertexIter::new(self)
    }

    //TODO iterate over references
    pub fn canonic_vertices(&self) -> impl Iterator<Item=Vertex> + '_ {
        (1..self.node_cnt()).map(|i| Vertex::forward(i))
    }

    pub fn link_cnt(&self) -> usize {
        self.all_links().count()
    }

    pub fn connector(&self, v: Vertex, w: Vertex) -> Option<Link> {
        //TODO rewrite via filter
        for l in self.outgoing_edges(v) {
            if l.end == w {
                return Some(l);
            }
        }
        None
    }

    pub fn v_str(&self, v: Vertex) -> String {
        format!("{}{}", self.node(v.node_id).name, Direction::str(v.direction))
    }

    pub fn gaf_str(&self, v: Vertex) -> String {
        format!("{}{}", Direction::gaf_str(v.direction), self.node(v.node_id).name)
    }

    pub fn l_str(&self, l: Link) -> String {
        format!("{}->{}", self.v_str(l.start), self.v_str(l.end))
    }

}

pub struct Path {
    v_storage: Vec<Vertex>,
    l_storage: Vec<Link>,
}

//FIXME rename, doesn't know about haplotype!
//Never empty! Use None instead
impl Path {

    pub fn new(init_v: Vertex) -> Path {
        Path {
            v_storage: vec![init_v],
            l_storage: Vec::new(),
        }
    }

    pub fn vertices(&self) -> &Vec<Vertex> {
        &self.v_storage
    }

    pub fn start(&self) -> Vertex {
        self.v_storage[0]
    }

    pub fn end(&self) -> Vertex {
        self.v_storage[self.v_storage.len() - 1]
    }

    pub fn len(&self) -> usize {
        self.v_storage.len()
    }

    pub fn links(&self) -> &Vec<Link> {
        &self.l_storage
    }

    pub fn reverse_complement(self) -> Path {
        //TODO optimize since consuming self
        Path {
            v_storage: self.v_storage.iter().rev().map(|v| v.rc()).collect(),
            l_storage: self.l_storage.iter().rev().map(|l| l.rc()).collect(),
        }
    }

    pub fn trim_to(&mut self, v: &Vertex) -> bool {
        if self.v_storage.contains(v) {
            while self.v_storage.last().unwrap() != v {
                self.v_storage.pop();
                self.l_storage.pop();
            }
            return true;
        }
        false
    }

    pub fn append(&mut self, l: Link) {
        assert!(self.v_storage.last().unwrap() == &l.start);
        //TODO disable expensive assert?
        assert!(!self.in_path(l.end.node_id));
        self.v_storage.push(l.end);
        self.l_storage.push(l);
    }

    pub fn in_path(&self, node_id: usize) -> bool {
        self.v_storage.iter().any(|v| v.node_id == node_id)
    }

    pub fn can_merge_in(&self, path: &Path) -> bool {
        assert!(self.v_storage.last() == path.v_storage.first());
        !path.v_storage.iter().skip(1).any(|v| self.in_path(v.node_id))
    }

    pub fn merge_in(&mut self, path: Path) {
        assert!(self.can_merge_in(&path));
        for l in path.l_storage {
            self.append(l);
        }
    }

    pub fn print(&self, g: &Graph) -> String {
        self.v_storage.iter().map(|&v| g.v_str(v)).collect::<Vec<String>>().join(",")
    }

    pub fn print_gaf(&self, g: &Graph) -> String {
        self.v_storage.iter().map(|&v| g.gaf_str(v)).collect::<Vec<String>>().join("")
    }

}
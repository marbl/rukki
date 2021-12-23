use std::str;
use std::collections::HashMap;

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
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
}

pub struct Node {
    //node size
    pub name: String,
    pub length: usize,
    pub coverage: f64,
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Vertex {
    //node id
    pub node_id: usize,
    //direction
    pub direction: Direction,
}

impl Vertex {
    fn rc(&self) -> Vertex {
        Vertex {
            node_id: self.node_id,
            direction: Direction::flip(self.direction),
        }
    }
}

//TODO separate 'links' and 'edges'
//links will have overlap size, CIGAR, etc
//edges will represent a Vertex pair
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Link {
    pub start: Vertex,
    pub end: Vertex,
    pub overlap: u32,
}

impl Link {
    fn rc(&self) -> Link {
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

impl Graph {

    fn new() -> Graph {
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

    fn node_iter(&self) -> std::slice::Iter<Node> {
        self.nodes.iter()
    }

    fn add_node(&mut self, node: Node) {
        //TODO rewrite without cloning with lifetimes
        self.name2ids.insert(node.name.clone(), self.nodes.len());
        self.nodes.push(node);
        self.incoming_links.push(Vec::new());
        self.outgoing_links.push(Vec::new());
    }

    fn add_link(&mut self, link: Link) {
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
        for node_id in 0..self.node_cnt() {
            let v = Vertex{node_id, direction: Direction::FORWARD};
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
            .nth(0)
    }

    fn parse_overlap(cigar: &str) -> u32 {
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
                let coverage: f64 = Self::parse_tag(&split[3..split.len()], "ll:f:")
                                        .unwrap_or(0.);
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

    fn get_vertex(&self, name: &str, direction: Direction) -> Vertex {
        let node_id = self.name2id(name);
        Vertex {node_id, direction}
    }

    fn rc(links: &Vec<Link>) -> Vec<Link> {
        links.into_iter().map(|x| x.rc()).collect()
    }

    pub fn node(&self, node_id: usize) -> &Node {
        &self.nodes[node_id]
    }

    //TODO switch to iterators when learn enough Rust :)
    fn outgoing_edges(&self, v: Vertex) -> Vec<Link> {
        match v.direction {
            Direction::FORWARD => self.outgoing_links[v.node_id].clone(),
            Direction::REVERSE => Self::rc(&self.incoming_links[v.node_id]),
        }
    }

    //TODO switch to iterators when learn enough Rust :)
    fn incoming_edges(&self, v: Vertex) -> Vec<Link> {
        match v.direction {
            Direction::FORWARD => self.incoming_links[v.node_id].clone(),
            Direction::REVERSE => Self::rc(&self.outgoing_links[v.node_id]),
        }
    }

    fn name2id(&self, name: &str) -> usize {
        match self.name2ids.get(name) {
            Some(&id) => id,
            None => panic!("Node {} is not in the graph", name),
        }
    }

    //FIXME figure out implicit lifetime
    pub fn all_links(&self) -> impl Iterator<Item=Link> + '_ {
        AllLinkIter::new(self)
    }

    pub fn all_nodes(&self) -> impl Iterator<Item=&Node> + '_ {
        self.nodes.iter()
    }

    pub fn link_cnt(&self) -> usize {
        self.all_links().count()
    }

    fn v_str(&self, v: Vertex) -> String {
        format!("{}{}", self.node(v.node_id).name, Direction::str(v.direction))
    }

    pub fn l_str(&self, l: Link) -> String {
        format!("{}->{}", self.v_str(l.start), self.v_str(l.end))
    }
}

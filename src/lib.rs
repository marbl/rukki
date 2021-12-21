use std::fs;
use std::env;
use std::str;
use std::error::Error;
use std::collections::HashMap;

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
enum Direction {
    FORWARD,
    REVERSE,
}

impl Direction {
    fn flip(d: Direction) -> Direction {
        match d {
            Direction::FORWARD => Direction::REVERSE,
            Direction::REVERSE => Direction::FORWARD,
        }
    }

    fn parse_char(c: char) -> Direction {
        match c {
            '+' => Direction::FORWARD,
            '-' => Direction::REVERSE,
            _ => panic!("Unknown direction {}", c),
        }
    }

    fn parse(s: &str) -> Direction {
        assert!(s.len() == 1, "Unknown direction {}", s);
        Direction::parse_char(s.chars().nth(0).unwrap())
    }
}

struct Node {
    //node size
    name: String,
    length: usize,
    coverage: f64,
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
struct Vertex {
    //node id
    node_id: usize,
    //direction
    direction: Direction,
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
#[derive(Copy, Clone, PartialEq, PartialOrd)]
struct Link {
    start: Vertex,
    end: Vertex,
    overlap: u32,
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
        Link::join_same(l1, l2) || Link::join_same(l1, &l2.rc())
    }
}

struct Graph {
    nodes: Vec<Node>,
    //incoming & outgoing links for every node
    incoming_links: Vec<Vec<Link>>,
    outgoing_links: Vec<Vec<Link>>,
    //TODO switch to &str and figure out how to work with lifetimes
    name2ids: HashMap<String, usize>,
}

impl Graph {

    fn new() -> Graph {
        let g = Graph {
            nodes: Vec::new(),
            incoming_links: Vec::new(),
            outgoing_links: Vec::new(),
            name2ids:  HashMap::new(),
        };

        if g.nodes.len() != g.incoming_links.len() || g.nodes.len() != g.outgoing_links.len() {
            panic!("Too bad")
        }
        g
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
        match link.end.direction {
            Direction::FORWARD => self.incoming_links[link.end.node_id].push(link),
            Direction::REVERSE => self.outgoing_links[link.end.node_id].push(link.rc()),
        };
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
    fn read(graph_str: &str) -> Graph {
        let mut g = Graph::new();

        for line in graph_str.lines() {
            if line.starts_with("S\t") {
                let split: Vec<&str> = line.split('\t').collect();
                //println!("Node line {:?}", split);
                let name = String::from(split[1]);
                let length = if split[2] != "*" {
                                 split[2].trim().len()
                             } else {
                                 Graph::parse_tag(&split[3..split.len()], "LN:i:")
                                     .expect("Neither sequence nor LN tag provided")
                             };
                let coverage: f64 = Graph::parse_tag(&split[3..split.len()], "ll:f:")
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
                let overlap = Graph::parse_overlap(split[5]);
                g.add_link(Link{start, end, overlap});
            }
        }
        g
    }

    fn get_vertex(&self, name: &str, direction: Direction) -> Vertex {
        let node_id = self.name2id(name);
        Vertex {node_id, direction}
    }

    fn rc(links: &Vec<Link>) -> Vec<Link> {
        links.into_iter().map(|x| x.rc()).collect()
    }

    fn node(&self, node_id: usize) -> &Node {
        &self.nodes[node_id]
    }

    //TODO switch to iterators when learn enough Rust :)
    fn outgoing_edges(&self, v: Vertex) -> Vec<Link> {
        match v.direction {
            Direction::FORWARD => self.outgoing_links[v.node_id].clone(),
            Direction::REVERSE => Graph::rc(&self.incoming_links[v.node_id]),
        }
    }

    //TODO switch to iterators when learn enough Rust :)
    fn incoming_edges(&self, v: Vertex) -> Vec<Link> {
        match v.direction {
            Direction::FORWARD => self.incoming_links[v.node_id].clone(),
            Direction::REVERSE => Graph::rc(&self.outgoing_links[v.node_id]),
        }
    }

    fn name2id(&self, name: &str) -> usize {
        match self.name2ids.get(name) {
            Some(&id) => id,
            None => panic!("Node {} is not in the graph", name),
        }
    }
}

pub struct Config {
    pub graph_fn: String,
}

impl Config {
    pub fn new(mut args: env::Args) -> Result<Config, &'static str> {
        //if args.len() < 2 {
        //    return Err("not enough arguments");
        //}
        args.next();

        let graph_fn = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get a graph path"),
        };

        Ok(Config { graph_fn })
    }
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let n = Node {
        name: String::from("Node 1"),
        length: 1,
        coverage: 0.,
    };
    let v = Vertex {
        node_id: 1,
        direction: Direction::FORWARD,
    };

    println!("Reading the graph from {:?}", &config.graph_fn);
    let g = Graph::read(&fs::read_to_string(&config.graph_fn)?);

    println!("Graph read successfully");

    Ok(())
}

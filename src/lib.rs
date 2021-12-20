use std::fs;
use std::env;
use std::error::Error;
use std::collections::HashMap;

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
enum Direction {
    FORWARD,
    REVERSE,
}

fn flip(d: Direction) -> Direction {
    match d {
        Direction::FORWARD => Direction::REVERSE,
        Direction::REVERSE => Direction::FORWARD,
    }
}

struct Node {
    //node size
    name: String,
    length: u64,
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
            direction: flip(self.direction),
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
    };
    let v = Vertex {
        node_id: 1,
        direction: Direction::FORWARD,
    };

    let g = Graph::new();

    println!("Hello, world! {:?}", v);

    let contents = fs::read_to_string(config.graph_fn)?;
    for line in contents.lines() {
        println!("Printing back {}", line)
    }
    Ok(())
}

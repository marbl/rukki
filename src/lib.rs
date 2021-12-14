use std::fs;
use std::error::Error;

#[derive(Copy, Clone, Debug)]
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

#[derive(Copy, Clone, Debug)]
struct Vertex {
    //node id
    node_id: usize,
    //direction
    direction: Direction,
    del: bool
}

impl Vertex {
    fn rc(&self) -> Vertex {
        Vertex {
            node_id: self.node_id,
            direction: flip(self.direction),
            del: self.del
        }
    }
}

#[derive(Copy, Clone)]
struct Link {
    start: Vertex,
    end: Vertex,
    overlap: u32,
    del: bool
}

impl Link {
    fn rc(&self) -> Link {
        Link {
            start: self.end.rc(),
            end: self.start.rc(),
            overlap: self.overlap,
            del: self.del,
        }
    }
}

struct Graph {
    nodes: Vec<Node>,
    //incoming & outgoing links for every node
    incoming_links: Vec<Vec<Link>>,
    outgoing_links: Vec<Vec<Link>>,
}

impl Graph {

    fn new() -> Graph {
        let g = Graph {
            nodes: Vec::new(),
            incoming_links: Vec::new(),
            outgoing_links: Vec::new(),
        };

        if g.nodes.len() != g.incoming_links.len() || g.nodes.len() != g.outgoing_links.len() {
            panic!("Too bad")
        }
        g
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
}

pub struct Config {
    pub graph_fn: String,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &str> {
        if args.len() < 2 {
            return Err("not enough arguments");
        }

        let graph_fn = args[1].clone();

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
        del: false
    };

    let g = Graph::new();

    println!("Hello, world! {:?}", v);

    let graph_str = fs::read_to_string(config.graph_fn)?;
    Ok(())
}

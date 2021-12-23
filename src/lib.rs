use std::fs;
use std::env;
use std::error::Error;

mod graph;

pub use graph::Graph;
pub use graph::Link;
pub use graph::Direction;

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
    println!("Reading the graph from {:?}", &config.graph_fn);
    let g = Graph::read(&fs::read_to_string(&config.graph_fn)?);

    println!("Graph read successfully");
    println!("Node count: {}", g.node_cnt());
    println!("Link count: {}", g.link_cnt());

    for l in g.all_links() {
        println!("here: {}", g.l_str(l));
    }

    Ok(())
}

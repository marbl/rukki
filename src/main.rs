//use std::io;
use std::env;
use std::process;
use graph_analysis::Config;

fn main() {
    let args: Vec<String> = env::args().collect();
    println!("Cmd arguments: {:?}", args);

    let config = Config::new(&args).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    println!("Processing graph in file {}", config.graph_fn);
    graph_analysis::run(config);
}

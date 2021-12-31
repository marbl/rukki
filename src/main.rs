//use std::io;
use std::env;
use std::process;
use graph_analysis::Config;
#[macro_use]
extern crate log;
use env_logger::{Builder, Target};

fn main() {
    //env_logger::init();
    let mut builder = Builder::from_default_env();
    builder.target(Target::Stdout);
    builder.init();
    info!("Starting up");
    //let args: Vec<String> = env::args().collect();
    println!("Cmd arguments: {:?}", env::args());

    let config = Config::new(env::args()).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    println!("Processing graph in file {}", config.graph_fn);
    match graph_analysis::run(config) {
        Ok(()) => println!("Success"),
        Err(e) => println!("Some error happened {:?}", e)
    }
}

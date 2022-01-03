//use std::io;
use std::env;
use graph_analysis::Config;
#[macro_use]
extern crate log;
use env_logger::{Env, Builder, Target};
use clap::Parser;

/// Assembly graph analysis
#[derive(Parser, Debug)]
#[clap(about, version, author)]
struct Args {
    /// GFA file
    //#[clap(short, long)]
    graph: String,

    /// Parental markers file
    #[clap(short, long)]
    parent_markers: String,

    /// Node annotation output file
    #[clap(short, long)]
    node_annotation: Option<String>,

    /// Marker-assisted extracted haplo-paths
    #[clap(short, long)]
    haplo_paths: Option<String>,
}

fn main() {
    //env_logger::init();
    let mut builder = Builder::from_env(Env::default().default_filter_or("info"));
    builder.target(Target::Stdout);
    builder.init();
    info!("Starting up");

    info!("Cmd arguments: {:?}", env::args());

    let args = Args::parse();

    //let config = Config::new(env::args()).unwrap_or_else(|err| {
    //    info!("Problem parsing arguments: {}", err);
    //    process::exit(1);
    //});
    let config = Config {
        graph_fn: args.graph,
        trio_markers_fn: args.parent_markers,
        init_node_annotation_fn: args.node_annotation,
        haplo_paths_fn: args.haplo_paths,
    };

    match graph_analysis::run(config) {
        Ok(()) => info!("Success"),
        Err(e) => info!("Some error happened {:?}", e)
    }
}

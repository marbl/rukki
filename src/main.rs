//use std::io;
use std::env;
#[macro_use]
extern crate log;
use env_logger::{Env, Builder, Target};
use clap::{AppSettings, Parser, Subcommand};

/// Assembly graph analysis suite
#[derive(Parser)]
#[clap(about, version, author)]
struct Args {
    /// GFA file
    //#[clap(short, long)]
    input_graph: String,

    //TODO add those
    ///// Sets a custom config file. Could have been an Option<T> with no default too
    //#[clap(short = "c", long = "config", default_value = "default.conf")]
    //config: String,

    // /// A level of verbosity, and can be used multiple times
    // #[clap(short = "v", long = "verbose", parse(from_occurrences))]
    // verbose: i32,

    #[clap(subcommand)]
    subcmd: Commands,

}

#[derive(Subcommand)]
enum Commands {
    /// Trio-marker based analysis
    #[clap(setting(AppSettings::ArgRequiredElseHelp))]
    Trio(TrioSettings),
    /// Primary-alt style analysis
    #[clap(setting(AppSettings::ArgRequiredElseHelp))]
    PriAlt,
}

//TODO use PathBuf?
#[derive(clap::Args)]
struct TrioSettings {
    /// Parental markers file
    #[clap(short, long)]
    parent_markers: String,

    /// Node annotation output file
    #[clap(short, long)]
    node_annotation: Option<String>,

    /// Marker-assisted extracted haplo-paths
    #[clap(long)]
    haplo_paths: Option<String>,

    /// Use GAF ([<>]<name1>)+ format for paths
    #[clap(short, long)]
    gaf_paths: bool,
}

fn main() {
    //env_logger::init();
    let mut builder = Builder::from_env(Env::default().default_filter_or("info"));
    builder.target(Target::Stdout);
    builder.init();
    info!("Starting up");

    info!("Cmd arguments: {:?}", env::args());

    let args = Args::parse();

    match &args.subcmd {
        Commands::Trio(settings) => {
            println!("Running trio marker analysis");

            match graph_analysis::run_trio_analysis(&args.input_graph, &settings.parent_markers,
                &settings.node_annotation, &settings.haplo_paths, settings.gaf_paths) {
                Ok(()) => info!("Success"),
                Err(e) => info!("Some error happened {:?}", e)
            }
        }

        Commands::PriAlt => {
            println!("Extracting primary/alt paths");
            match graph_analysis::run_primary_alt_analysis(&args.input_graph) {
                Ok(()) => info!("Success"),
                Err(e) => info!("Some error happened {:?}", e)
            }
        }
    }

}

//use std::io;
#[macro_use]
extern crate log;
use env_logger::{Env, Builder, Target};
use clap::{AppSettings, Parser, Subcommand};

/// Assembly graph analysis suite
#[derive(Parser)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(about, version, author)]
struct Args {
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
    Trio(graph_analysis::TrioSettings),
    // /// Primary-alt style analysis
    // #[clap(setting(AppSettings::ArgRequiredElseHelp))]
    // PriAlt(PriAltSettings),
}

//TODO use PathBuf?
#[derive(clap::Args)]
struct PriAltSettings {
    /// GFA file
    #[clap(short, long)]
    graph: String,

    /// colors output file
    #[clap(short, long)]
    assign: Option<String>,

    /// Extracted primary and alt paths
    #[clap(long)]
    paths: Option<String>,

    /// Use GAF ([<>]utg_name)+ format for paths
    #[clap(short, long)]
    gaf_format: bool,

}

fn main() {
    //env_logger::init();
    let mut builder = Builder::from_env(Env::default().default_filter_or("info"));
    builder.target(Target::Stdout);
    builder.init();
    //info!("Starting up");

    //info!("Cmd arguments: {:?}", env::args());

    let args = Args::parse();

    match &args.subcmd {
        Commands::Trio(settings) => {
            info!("Running trio marker analysis");

            match graph_analysis::run_trio_analysis(&settings) {
                Ok(()) => info!("Success"),
                Err(e) => info!("Some error happened {:?}", e)
            }
        }

        //Commands::PriAlt(settings) => {
        //    info!("Extracting primary/alt paths");
        //    match graph_analysis::run_primary_alt_analysis(&settings.graph,
        //                                &settings.assign, &settings.paths,
        //                                settings.gaf_paths) {
        //        Ok(()) => info!("Success"),
        //        Err(e) => info!("Some error happened {:?}", e)
        //    }
        //}
    }

}

//use std::io;
#[macro_use]
extern crate log;
use clap::{Parser, Subcommand};
use env_logger::{Builder, Env, Target};

#[derive(Parser, Debug)]
#[command(name = "rukki", author = "Sergey Nurk", about = "extraction of paths from assembly graphs", long_about=None)]
struct Args {
    #[clap(subcommand)]
    subcmd: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Trio-marker based analysis
    Trio(rukki::TrioSettings),
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
            settings.validate();

            match rukki::run_trio_analysis(settings) {
                Ok(()) => info!("Success"),
                Err(e) => info!("Some error happened {:?}", e),
            }
        }
    }
}

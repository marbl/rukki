use std::fs;
use std::env;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use log::info;

mod graph;
mod trio;
mod trio_walk;

pub use graph::Graph;
pub use graph::Vertex;
pub use graph::Link;
pub use graph::Direction;

pub struct Config {
    pub graph_fn: String,
    pub trio_markers_fn: String,
    pub out_fn: String,
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

        let trio_markers_fn = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get trio marker path"),
        };

        let out_fn = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get an output path"),
        };

        Ok(Config { graph_fn, trio_markers_fn, out_fn })
    }
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    info!("Reading the graph from {:?}", &config.graph_fn);
    let g = Graph::read(&fs::read_to_string(&config.graph_fn)?);

    info!("Graph read successfully");
    info!("Node count: {}", g.node_cnt());
    info!("Link count: {}", g.link_cnt());

    //for n in g.all_nodes() {
    //    println!("Node: {} length: {} cov: {}", n.name, n.length, n.coverage);
    //}

    //for l in g.all_links() {
    //    println!("Link: {}", g.l_str(l));
    //}

    let trio_infos = trio::read_trio(&fs::read_to_string(&config.trio_markers_fn)?);
    let parental_groups = trio::assign_parental_groups(&g, &trio_infos);

    let mut output = File::create(config.out_fn)?;
    //write!(output, "{}", g.as_gfa())?;

    //for (node_id, n) in g.all_nodes().enumerate() {
    //    assert!(g.name2id(&n.name) == node_id);
    //    println!("Node: {} length: {} cov: {}", n.name, n.length, n.coverage);
    //    match parental_groups.get(node_id) {
    //        None => println!("No assignment"),
    //        Some(assign) => println!("Assinged. Info: {:?}", assign),
    //    }
    //}

    writeln!(output, "node\tmat:pat\tassignment\tcolor")?;
    for (node_id, n) in g.all_nodes().enumerate() {
        assert!(g.name2id(&n.name) == node_id);
        if let Some(assign) = parental_groups.get(node_id) {
            let color = match assign.group {
                trio::TrioGroup::PATERNAL => "#8888FF",
                trio::TrioGroup::MATERNAL => "#FF8888",
                trio::TrioGroup::ISSUE => "#fbb117",
                trio::TrioGroup::HOMOZYGOUS => "#c5d165",
            };
            writeln!(output, "{}\t{}\t{:?}\t{}", n.name, assign.info, assign.group, color)?;
        }
    }

    let init_node_len_thr = 500_000;
    let mut path_searcher = trio_walk::HaploPathSearcher::new(&g,
        &parental_groups, init_node_len_thr);

    for (path, group) in path_searcher.find_all() {
        info!("Identified {:?} path: {}", group, path.print(&g));
    }

    Ok(())
}
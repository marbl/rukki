use std::fs;
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
    pub init_node_annotation_fn: Option<String>,
    pub haplo_paths_fn: Option<String>,
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    info!("Reading graph from {:?}", &config.graph_fn);
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
    //write!(output, "{}", g.as_gfa())?;

    info!("Reading trio marker information from {}", config.trio_markers_fn);
    let trio_infos = trio::read_trio(&fs::read_to_string(&config.trio_markers_fn)?);

    info!("Assigning initial parental groups to the nodes");
    let parental_groups = trio::assign_parental_groups(&g, &trio_infos);

    if let Some(output) = config.init_node_annotation_fn {
        info!("Writing initial node annotation to {}", output);
        let mut output = File::create(output)?;

        writeln!(output, "node\tlength\tmat:pat\tassignment\tcolor")?;
        for (node_id, n) in g.all_nodes().enumerate() {
            assert!(g.name2id(&n.name) == node_id);
            if let Some(assign) = parental_groups.get(node_id) {
                let color = match assign.group {
                    trio::TrioGroup::PATERNAL => "#8888FF",
                    trio::TrioGroup::MATERNAL => "#FF8888",
                    trio::TrioGroup::ISSUE => "#fbb117",
                    trio::TrioGroup::HOMOZYGOUS => "#c5d165",
                };
                writeln!(output, "{}\t{}\t{}\t{:?}\t{}", n.name, n.length, assign.info
                                                       , assign.group, color)?;
            }
        }
    }

    if let Some(output) = config.haplo_paths_fn {
        info!("Searching for haplo-paths, output in {}", output);
        let mut output = File::create(output)?;

        writeln!(output, "name\tpath\tassignment\tinit_node")?;
        let init_node_len_thr = 500_000;
        let mut path_searcher = trio_walk::HaploPathSearcher::new(&g,
            &parental_groups, init_node_len_thr);

        for (path, group) in path_searcher.find_all() {
            assert!(path.vertices().contains(&Vertex{node_id: path.initial_node(),
                                                     direction: Direction::FORWARD}));
            //info!("Identified {:?} path: {}", group, path.print(&g));
            writeln!(output, "path_from_{}\t{}\t{:?}\t{}",
                g.node(path.initial_node()).name,
                path.print(&g),
                group,
                path.initial_node())?;
        }

        let used = path_searcher.used();

        for (node_id, n) in g.all_nodes().enumerate() {
            if !used.contains_key(&node_id) {
                let group_str = parental_groups.group(node_id)
                                    .map_or(String::from("NA"), |x| format!("{:?}", x));

                //println!("Unused node: {} length: {} group: {}", n.name, n.length, group_str);
                writeln!(output, "unused_{}\t{}\t{}\t{}",
                    n.name,
                    &g.v_str(Vertex::forward(node_id)),
                    group_str,
                    node_id)?;
            }
        }

    }

    info!("All done");
    Ok(())
}
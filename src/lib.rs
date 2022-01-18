use std::fs;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use log::info;

mod graph;
//ASK why tests don't compile without the pub
pub mod graph_algos;
mod trio;
mod trio_walk;

pub use graph::Graph;
pub use graph::Vertex;
pub use graph::Link;
pub use graph::Direction;

pub fn run_trio_analysis(graph_fn: &String, trio_markers_fn: &String,
    init_node_annotation_fn: &Option<String>, haplo_paths_fn: &Option<String>,
    gaf_paths: bool) -> Result<(), Box<dyn Error>> {
    info!("Reading graph from {:?}", graph_fn);
    let g = Graph::read(&fs::read_to_string(graph_fn)?);

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

    info!("Reading trio marker information from {}", trio_markers_fn);
    let trio_infos = trio::read_trio(&fs::read_to_string(trio_markers_fn)?);

    info!("Assigning initial parental groups to the nodes");
    let parental_groups = trio::assign_parental_groups(&g, &trio_infos);
    info!("Detecting homozygous nodes");
    //TODO parameterize
    let parental_groups = trio_walk::assign_homozygous(&g, parental_groups, 100_000);

    if let Some(output) = init_node_annotation_fn {
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

    if let Some(output) = haplo_paths_fn {
        info!("Searching for haplo-paths, output in {}", output);
        let mut output = File::create(output)?;

        writeln!(output, "name\tpath\tassignment\tinit_node")?;
        let init_node_len_thr = 500_000;
        let mut path_searcher = trio_walk::HaploSearcher::new(&g,
            &parental_groups, init_node_len_thr);

        for (path, node_id, group) in path_searcher.find_all() {
            assert!(path.vertices().contains(&Vertex::forward(node_id)));
            //info!("Identified {:?} path: {}", group, path.print(&g));
            writeln!(output, "path_from_{}\t{}\t{:?}\t{}",
                g.node(node_id).name,
                if gaf_paths {path.print_gaf(&g)}
                               else {path.print(&g)},
                group,
                g.node(node_id).name)?;
        }

        let used = path_searcher.used();

        for (node_id, n) in g.all_nodes().enumerate() {
            if !used.contains_key(&node_id) {
                let group_str = parental_groups.group(node_id)
                                    .map_or(String::from("NA"), |x| format!("{:?}", x));

                //println!("Unused node: {} length: {} group: {}", n.name, n.length, group_str);
                writeln!(output, "unused_{}_len_{}\t{}\t{}\t{}",
                    n.name,
                    n.length,
                    if gaf_paths {g.gaf_str(Vertex::forward(node_id))}
                                else {g.v_str(Vertex::forward(node_id))},
                    group_str,
                    node_id)?;
            }
            //FIXME how many times should we report HOMOZYGOUS node?!
            //What if it has never been used? Are we confident enough?
        }

    }

    info!("All done");
    Ok(())
}
use std::fs;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use log::info;
use std::collections::HashSet;

//tests don't compile without the pub
//FIXME what to do?
pub mod graph;
pub mod graph_algos;
pub mod trio;
pub mod trio_walk;
pub mod pseudo_hap;

pub use graph::Graph;
pub use graph::Vertex;
pub use graph::Path;
pub use graph::Link;
pub use graph::Direction;

fn read_graph(graph_fn: &str) -> Result<Graph, Box<dyn Error>>  {
    info!("Reading graph from {:?}", graph_fn);
    let g = Graph::read_sanitize(&fs::read_to_string(graph_fn)?);

    info!("Graph read successfully");
    info!("Node count: {}", g.node_cnt());
    info!("Link count: {}", g.link_cnt());
    Ok(g)
}

fn output_coloring(g: &Graph,
                   assignments: &trio::AssignmentStorage,
                   file_name: &str)
                   -> Result<(), std::io::Error> {
    let mut output = File::create(file_name)?;
    writeln!(output, "node\tassignment\tlength\tinfo\tcolor")?;
    for (node_id, n) in g.all_nodes().enumerate() {
        assert!(g.name2id(&n.name) == node_id);
        if let Some(assign) = assignments.get(node_id) {
            let color = match assign.group {
                trio::TrioGroup::PATERNAL => "#8888FF",
                trio::TrioGroup::MATERNAL => "#FF8888",
                trio::TrioGroup::ISSUE => "#FFDE24",
                trio::TrioGroup::HOMOZYGOUS => "#7900D6",
            };
            writeln!(output, "{}\t{:?}\t{}\t{}\t{}", n.name
                                                   , assign.group
                                                   , n.length
                                                   , assign.info
                                                   , color)?;
        }
    }
    Ok(())
}

pub fn augment_by_path_search<'a>(g: &'a Graph,
    assignments: &'a trio::AssignmentStorage<'a>,
    init_node_len_thr: usize) -> trio::AssignmentStorage<'a> {

    let mut assigning_path_searcher = trio_walk::HaploSearcher::new_assigning(&g,
        &assignments, init_node_len_thr);

    assigning_path_searcher.find_all();
    let tentative_assignments = assigning_path_searcher.used();

    let mut init_assign = assignments.clone();
    for node_id in tentative_assignments.assigned() {
        let tentative_group = tentative_assignments.group(node_id).unwrap();
        //any mixed assignment has chance to be erroneous due to graph issues
        if !tentative_group.is_definite() {
            continue;
        }
        match init_assign.group(node_id) {
            None => {
                info!("Assigning tentative group {:?} to node {}", tentative_group, g.name(node_id));
                init_assign.assign(node_id, trio::Assignment::<trio::TrioGroup> {
                    group: tentative_group,
                    confidence: trio::Confidence::MODERATE,
                    info: String::from("PreliminaryLaunch"),
                });},
            Some(init_group) => assert!(init_group == tentative_group
                    || init_group == trio::TrioGroup::HOMOZYGOUS),
        }
    }
    init_assign
}

pub fn run_trio_analysis(graph_fn: &str, trio_markers_fn: &str,
    init_node_annotation_fn: &Option<String>,
    final_node_annotation_fn: &Option<String>,
    haplo_paths_fn: &Option<String>,
    gaf_paths: bool, low_cnt_thr: usize, ratio_thr: f32,
    low_inv_density: usize) -> Result<(), Box<dyn Error>> {
    let g = read_graph(graph_fn)?;

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
    let init_assign = trio::assign_parental_groups(&g, &trio_infos,
        low_cnt_thr, ratio_thr, low_inv_density);
    //TODO parameterize
    let init_assign = trio::assign_homozygous(&g, init_assign, 100_000);
    //let init_assign = trio::assign_small_bubbles(&g, init_assign, 100_000);

    if let Some(output) = init_node_annotation_fn {
        info!("Writing initial node annotation to {}", output);
        output_coloring(&g, &init_assign, output)?;
    }

    let init_node_len_thr = 500_000;
    //FIXME parameterize
    let augment_assign = augment_by_path_search(&g,
        &init_assign,
        init_node_len_thr);

    let mut path_searcher = trio_walk::HaploSearcher::new(&g,
        &augment_assign, init_node_len_thr);

    let haplo_paths = path_searcher.find_all();

    let mut final_assign = path_searcher.used().clone();

    for (node_id, _) in g.all_nodes().enumerate() {
        if !final_assign.contains(node_id) && init_assign.contains(node_id) {
            //FIXME how many times should we report HOMOZYGOUS node?!
            //What if it has never been used? Are we confident enough?
            final_assign.assign(node_id, init_assign.get(node_id).unwrap().clone())
        }
    }

    if let Some(output) = haplo_paths_fn {
        info!("Searching for haplo-paths, output in {}", output);
        let mut output = File::create(output)?;

        writeln!(output, "name\tpath\tassignment\tinit_node")?;
        for (path, node_id, group) in haplo_paths {
            assert!(path.vertices().contains(&Vertex::forward(node_id)));
            //info!("Identified {:?} path: {}", group, path.print(&g));
            writeln!(output, "path_from_{}\t{}\t{:?}\t{}",
                g.node(node_id).name,
                path.print_format(&g, gaf_paths),
                group,
                g.node(node_id).name)?;
        }

        for (node_id, n) in g.all_nodes().enumerate() {
            if !path_searcher.used().contains(node_id) {
                let group_str = init_assign
                    .group(node_id)
                    .map_or(String::from("NA"), |x| format!("{:?}", x));

                //println!("Unused node: {} length: {} group: {}", n.name, n.length, group_str);
                writeln!(output, "unused_{}_len_{}\t{}\t{}\t{}",
                    n.name,
                    n.length,
                    Path::new(Vertex::forward(node_id)).print_format(&g, gaf_paths),
                    group_str,
                    node_id)?;
            }
            //FIXME how many times should we report HOMOZYGOUS node?!
            //What if it has never been used? Are we confident enough?
        }
    }

    if let Some(output) = final_node_annotation_fn {
        info!("Writing final node annotation to {}", output);
        output_coloring(&g, &final_assign, output)?;
    }

    info!("All done");
    Ok(())
}

pub fn run_primary_alt_analysis(graph_fn: &str,
                                colors_fn: &Option<String>,
                                paths_fn: &Option<String>,
                                gaf_paths: bool) -> Result<(), Box<dyn Error>> {
    let g = read_graph(graph_fn)?;
    let unique_block_len = 500_000;
    let linear_blocks = pseudo_hap::pseudo_hap_decompose(&g, unique_block_len);

    if let Some(output) = colors_fn {
        info!("Writing node colors to {}", output);
        let mut output = File::create(output)?;

        let mut primary_nodes = HashSet::new();
        let mut alt_nodes = HashSet::new();
        let mut boundary_nodes = HashSet::new();

        for block in &linear_blocks {
            let p = block.instance_path();
            primary_nodes.extend(p.vertices()
                                 .iter().map(|&v| v.node_id));
            alt_nodes.extend(block.known_alt_nodes().iter().copied());
            boundary_nodes.extend([p.start().node_id, p.end().node_id]);
        }

        writeln!(output, "node\tlength\tassignment\tcolor")?;
        for (node_id, n) in g.all_nodes().enumerate() {
            assert!(g.name2id(&n.name) == node_id);
            let mut color = "#808080";
            let mut assign = "NA";
            if boundary_nodes.contains(&node_id) {
                assert!(!alt_nodes.contains(&node_id));
                color = "#fbb117";
                assign = "PRIMARY_BOUNDARY";
            } else if primary_nodes.contains(&node_id) {
                assert!(!alt_nodes.contains(&node_id));
                color = "#8888FF";
                assign = "PRIMARY";
            } else if alt_nodes.contains(&node_id) {
                color = "#FF8888";
                assign = "ALT";
            }
            writeln!(output, "{}\t{}\t{}\t{}", n.name, n.length, assign, color)?;
        }
    }

    let used : HashSet<usize> = linear_blocks.iter()
                                    .flat_map(|b| b.all_nodes())
                                    .collect();

    if let Some(output) = paths_fn {
        info!("Outputting paths in {}", output);
        let mut output = File::create(output)?;

        writeln!(output, "name\tlen\tpath\tassignment")?;

        for (block_id, block) in linear_blocks.into_iter().enumerate() {
            writeln!(output, "primary_{}\t{}\t{}\tPRIMARY",
                block_id,
                block.instance_path().total_length(&g),
                block.instance_path().print_format(&g, gaf_paths))?;
            for (alt_id, &known_alt) in block.known_alt_nodes().iter().enumerate() {
                writeln!(output, "alt_{}_{}\t{}\t{}\tALT",
                    block_id,
                    alt_id,
                    g.node(known_alt).length,
                    Path::new(Vertex::forward(known_alt)).print_format(&g, gaf_paths))?;
            }
        }

        for (node_id, n) in g.all_nodes().enumerate() {
            if !used.contains(&node_id) {
                writeln!(output, "unused_{}\t{}\t{}\tNA",
                    n.name,
                    n.length,
                    Path::new(Vertex::forward(node_id)).print_format(&g, gaf_paths))?;
            }
        }
    }

    info!("All done");
    Ok(())
}

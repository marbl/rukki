use std::fs;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use log::{info,warn,debug};
use trio_walk::HaploSearchSettings;
use std::collections::HashSet;

//tests don't compile without the pub
//FIXME what to do?
pub mod graph;
pub mod graph_algos;
pub mod trio;
pub mod trio_walk;
pub mod pseudo_hap;

pub use graph::*;

use crate::trio::TrioGroup;

//TODO use PathBuf
#[derive(clap::Args)]
pub struct TrioSettings {
    /// GFA file
    #[clap(short, long)]
    pub graph: String,

    /// Parental markers file
    #[clap(short, long)]
    pub markers: String,

    /// Marker-based annotation output file
    #[clap(long)]
    pub init_assign: Option<String>,

    /// Refined annotation output file
    #[clap(long)]
    pub refined_assign: Option<String>,

    /// Final annotation output file
    #[clap(long)]
    pub final_assign: Option<String>,

    /// Marker-assisted extracted haplo-paths
    #[clap(long, short)]
    pub paths: Option<String>,

    /// Use GAF ([<>]<name1>)+ format for paths
    #[clap(long)]
    pub gaf_format: bool,

    /// Minimal number of parent-specific markers required for assigning parental group to a node
    #[clap(long, default_value_t = 10)]
    pub marker_cnt: usize,

    /// Require at least (node_length / <value>) markers within the node for parental group assignment
    #[clap(long, default_value_t = 10_000)]
    pub marker_sparsity: usize,

    /// Sets minimal marker excess for assigning a parental group to <value>:1
    #[clap(long, default_value_t = 5.0)]
    pub marker_ratio: f64,

    /// Longer nodes are unlikely to be spurious and likely to be reliably assigned based on markers (used in HOMOZYGOUS node labeling)
    #[clap(long, default_value_t = 200_000)]
    pub trusted_len: usize,

    /// Nodes with coverage below <coeff> * <weighted mean coverage of 'solid' nodes> can not be 'reclassified' as homozygous (negative turns off reclassification, 0. disables coverage check)
    #[clap(long, default_value_t = 1.5)]
    pub suspect_homozygous_cov_coeff: f64,

    /// Longer nodes can not be classified as homozygous
    #[clap(long, default_value_t = 1_000_000)]
    pub homozygous_max_len: usize,

    //TODO maybe check that it is > trusted_len
    /// Longer nodes are unlikely to represent repeats, polymorphic variants, etc (used to seed and guide the path search)
    #[clap(long, default_value_t = 500_000)]
    pub solid_len: usize,

    /// Minimal node length for assigning ISSUE label
    #[clap(long, default_value_t = 50_000)]
    pub issue_len: usize,

    /// Minimal number of markers for assigning ISSUE label (by default == marker_cnt, will typically be set to a value >= marker_cnt)
    #[clap(long)]
    pub issue_cnt: Option<usize>,

    /// Require at least (node_length / <value>) markers for assigning ISSUE label (by default == marker_sparsity, will typically be set to a value >= marker_sparsity)
    #[clap(long)]
    pub issue_sparsity: Option<usize>,

    /// Require primary marker excess BELOW <value>:1 for assigning ISSUE label. Must be <= marker_ratio (by default == marker_ratio)
    #[clap(long)]
    pub issue_ratio: Option<f64>,

    /// Try to fill in small ambiguous bubbles
    #[clap(long)]
    pub try_fill_bubbles: bool,

    /// Do not fill bubble if source or sink is non-solid, non-homozygous and has coverage above <coeff> * <weighted mean coverage of 'solid' nodes>
    /// (0. fails coverage check, negative makes it to always pass)
    #[clap(long, default_value_t = 1.5)]
    pub max_unique_cov_coeff: f64,

    /// Bubbles including a longer alternative sequence will not be filled
    #[clap(long, default_value_t = 50_000)]
    fillable_bubble_len: usize,

    /// Bubbles with bigger difference between alternatives' lengths will not be filled
    #[clap(long, default_value_t = 200)]
    fillable_bubble_diff: usize,

    /// Minimal introducible gap size (number of Ns reported). If the gap size estimate is smaller it will be artificially increased to this value.
    #[clap(long, default_value_t = 1000)]
    min_gap_size: usize,

    /// Default gap size, which will be output in cases where reasonable estimate is not possible or (more likely) hasn't been implemented yet.
    #[clap(long, default_value_t = 5000)]
    default_gap_size: usize,

}

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

pub fn augment_by_path_search(g: &Graph,
    assignments: trio::AssignmentStorage,
    settings: HaploSearchSettings) -> trio::AssignmentStorage {
    info!("Augmenting node annotation by path search. Round 1.");
    let assignments = augment_by_path_search_round(&g,
        assignments,
        settings);
    info!("Augmenting node annotation by path search. Round 2.");
    augment_by_path_search_round(&g,
        assignments,
        settings)
}

fn augment_by_path_search_round(g: &Graph,
    assignments: trio::AssignmentStorage,
    settings: HaploSearchSettings) -> trio::AssignmentStorage {

    let mut path_searcher = trio_walk::HaploSearcher::new(&g,
        &assignments, settings.assigning_stage_adjusted());

    path_searcher.find_all();
    let node_usage = path_searcher.take_used();
    augment_assignments(g, assignments, &node_usage, true)
}

fn augment_assignments(g: &Graph,
    mut assignments: trio::AssignmentStorage,
    extra_assignments: &trio::AssignmentStorage,
    exclude_homozygous: bool) -> trio::AssignmentStorage {
    for node_id in extra_assignments.assigned() {
        let tentative_group = extra_assignments.group(node_id).unwrap();
        assert!(tentative_group != TrioGroup::ISSUE);
        //any mixed assignment has chance to be erroneous due to graph issues
        if exclude_homozygous && !tentative_group.is_definite() {
            continue;
        }
        match assignments.group(node_id) {
            None => {
                debug!("Assigning tentative group {:?} to node {}", tentative_group, g.name(node_id));
                assignments.assign(node_id, tentative_group, String::from("PathSearch"));
            },
            Some(init_group) => assert!(init_group == tentative_group
                    || init_group == trio::TrioGroup::HOMOZYGOUS),
        }
    }
    assignments
}

fn weighted_mean_solid_cov(g: &Graph, solid_len_thr: usize) -> f64 {
    let mut total_len = 0;
    let mut total_cov = 0.;
    for n in g.all_nodes() {
        if n.length >= solid_len_thr {
            total_len += n.length;
            total_cov += n.coverage * (n.length as f64);
        }
    }
    total_cov / total_len as f64
}

pub fn run_trio_analysis(settings: &TrioSettings) -> Result<(), Box<dyn Error>> {
    let g = read_graph(&settings.graph)?;

    //for n in g.all_nodes() {
    //    println!("Node: {} length: {} cov: {}", n.name, n.length, n.coverage);
    //}
    //for l in g.all_links() {
    //    println!("Link: {}", g.l_str(l));
    //}
    //write!(output, "{}", g.as_gfa())?;

    info!("Reading trio marker information from {}", &settings.markers);
    let trio_infos = trio::read_trio(&fs::read_to_string(&settings.markers)?);

    info!("Assigning initial parental groups to the nodes");
    let assignments = trio::assign_parental_groups(&g, &trio_infos,
        settings.marker_cnt, settings.marker_sparsity, settings.marker_ratio,
        settings.issue_len,
        settings.issue_cnt.unwrap_or(settings.marker_cnt),
    settings.issue_sparsity.unwrap_or(settings.marker_sparsity),
        settings.issue_ratio.unwrap_or(settings.marker_ratio),
    );

    if let Some(output) = &settings.init_assign {
        info!("Writing initial node annotation to {}", output);
        output_coloring(&g, &assignments, output)?;
    }

    let solid_cov_est = weighted_mean_solid_cov(&g, settings.solid_len);
    let suspect_homozygous_cov = if settings.suspect_homozygous_cov_coeff <= 0. {
        settings.suspect_homozygous_cov_coeff
    } else {
        if solid_cov_est == 0. {
            warn!("Looks like the graph didn't have coverage information, which we were hoping to use. Consider providing it or changing --suspect-homozygous-cov-coeff");
        }
        settings.suspect_homozygous_cov_coeff * solid_cov_est
    };

    let assignments = trio::assign_homozygous(&g, assignments,
        settings.trusted_len, suspect_homozygous_cov,
        settings.homozygous_max_len);

    let mut search_settings = HaploSearchSettings {
        solid_len: settings.solid_len,
        trusted_len: settings.trusted_len,
        fillable_bubble_len: settings.fillable_bubble_len,
        fillable_bubble_diff: settings.fillable_bubble_diff,
        min_gap_size: settings.min_gap_size as i64,
        default_gap_size: settings.default_gap_size as i64,
        ..HaploSearchSettings::default()
    };

    if settings.try_fill_bubbles {
        search_settings.ambig_filling_level = 2;
        info!("Will try filling small bubbles");
        //assert!(settings.max_unique_cov_coeff >= 0.);
        if settings.max_unique_cov_coeff < 0. {
            //leaving default
            search_settings.max_unique_cov = f64::MAX;
            info!("Negative '--max-unique-cov-coeff' provided. All nodes will be considered unique for purposes of bubble filling");
        }
        if settings.max_unique_cov_coeff > 0. && solid_cov_est == 0. {
            warn!("Looks like the graph didn't have coverage information, which we were hoping to use. Consider providing it or changing --max-unique-cov-coeff");
        }
        search_settings.max_unique_cov = settings.max_unique_cov_coeff * solid_cov_est;
        info!("Maximal 'unique' coverage for bubble filling set to {}", search_settings.max_unique_cov);
        if search_settings.max_unique_cov == 0. {
            info!("Will only fill bubbles between solid or homozygous nodes");
        }
    }

    let assignments = augment_by_path_search(&g,
        assignments,
        search_settings);

    if let Some(output) = &settings.refined_assign {
        info!("Writing refined node annotation to {}", output);
        output_coloring(&g, &assignments, output)?;
    }
    let mut path_searcher = search_settings.build_searcher(&g, &assignments);

    let haplo_paths = path_searcher.find_all();
    let node_usage = path_searcher.take_used();

    let assignments = augment_assignments(&g, assignments,
        &node_usage, false);

    if let Some(output) = &settings.final_assign {
        info!("Writing final node annotation to {}", output);
        output_coloring(&g, &assignments, output)?;
    }

    let short_group_str = |o_g: Option<TrioGroup>| {
        match o_g {
            Some(TrioGroup::MATERNAL) => "mat",
            Some(TrioGroup::PATERNAL) => "pat",
            _ => "na",
        }
    };

    if let Some(output) = &settings.paths {
        info!("Outputting haplo-paths to {}", output);
        let mut output = File::create(output)?;

        writeln!(output, "name\tpath\tassignment")?;
        for (path, node_id, group) in haplo_paths {
            assert!(path.vertices().contains(&Vertex::forward(node_id)));
            //info!("Identified {:?} path: {}", group, path.print(&g));
            writeln!(output, "{}_from_{}\t{}\t{:?}",
                short_group_str(Some(group)),
                g.node(node_id).name,
                path.print_format(&g, settings.gaf_format),
                group)?;
        }

        let mut write_node = |n: &Node, group: Option<TrioGroup>| {
            writeln!(output, "{}_unused_{}\t{}\t{}",
                short_group_str(group),
                n.name,
                Direction::format_node(&n.name, Direction::FORWARD, settings.gaf_format),
                group.map_or(String::from("NA"), |x| format!("{:?}", x)))
        };

        for (node_id, n) in g.all_nodes().enumerate() {
            let haplopath_assign = node_usage.group(node_id);
            match assignments.group(node_id) {
                None | Some(TrioGroup::ISSUE) => {
                    assert!(!node_usage.contains(node_id));
                    debug!("Node: {} length: {} not assigned to any haplotype (adding trivial NA path)",
                        n.name, n.length);
                    write_node(g.node(node_id), None)?;
                }
                Some(assign) => {
                    if TrioGroup::compatible(assign, TrioGroup::MATERNAL)
                        //not present in haplopaths paths or incompatible
                        && haplopath_assign.map_or(true,
                            |x| TrioGroup::incompatible(x, TrioGroup::MATERNAL)) {
                        debug!("Node: {} length: {} not present in MATERNAL haplo-paths (adding trivial MATERNAL path)",
                            n.name, n.length);
                        write_node(g.node(node_id), Some(TrioGroup::MATERNAL))?;
                    }
                    if TrioGroup::compatible(assign, TrioGroup::PATERNAL)
                        //not present in haplopaths paths or incompatible
                        && haplopath_assign.map_or(true,
                            |x| TrioGroup::incompatible(x, TrioGroup::PATERNAL)) {
                        debug!("Node: {} length: {} not present in PATERNAL haplo-paths (adding trivial PATERNAL path)",
                            n.name, n.length);
                        write_node(g.node(node_id), Some(TrioGroup::PATERNAL))?;
                    }
                }
            }
        }
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

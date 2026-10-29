use log::{debug, info, warn};
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::{
    collections::{HashSet, VecDeque},
    path::PathBuf,
};
use trio_walk::HaploSearchSettings;

//tests don't compile without the pub
//FIXME what to do?
pub mod graph;
pub mod graph_algos;
pub mod pseudo_hap;
pub mod trio;
pub mod trio_walk;

pub use graph::*;

use crate::trio::{
    assign_short_node_tangles, GroupAssignmentSettings, TangleAssignmentSettings, TrioGroup,
};
use crate::trio_walk::HaploSearcher;

/// Settings for evidence-free dual assembly using the pseudo-hap bubble
/// decomposition as a deterministic source of complementary choices.
#[derive(clap::Args, Debug)]
pub struct DualSettings {
    /// GFA file with per-node coverage
    #[clap(short, long)]
    graph: PathBuf,

    /// Node coloring output
    #[clap(long)]
    colors: Option<PathBuf>,

    /// Extracted pseudohaplotype paths
    #[clap(long, short)]
    paths: Option<PathBuf>,

    /// Strict homologous node pairs selected from MashMap self-alignments
    #[clap(long)]
    homology: Option<PathBuf>,

    /// Use GAF ([<>]<name1>)+ format for paths
    #[clap(long)]
    gaf_format: bool,

    /// Comma separated names for the two outputs
    #[clap(long, default_value_t = String::from("haplotype1,haplotype2"))]
    hap_names: String,

    /// Seed used to swap primary and alternate colors between blocks
    #[clap(long, default_value_t = 1)]
    seed: u64,

    /// Treat graph-wide coverage as haploid depth when no bubble-based
    /// coverage estimate is available
    #[clap(long)]
    haploid: bool,

    /// Minimum length of a unique block used by pseudo-hap decomposition
    #[clap(long, default_value_t = 500_000)]
    unique_block_len: usize,

    /// Components below this multiple of haploid coverage are emitted once;
    /// components near diploid coverage are shared by both outputs.
    #[clap(long, default_value_t = 1.5)]
    haploid_component_cov_coeff: f64,
}

impl DualSettings {
    pub fn validate(&self) {
        assert!(self.unique_block_len > 0);
        assert!(self.haploid_component_cov_coeff > 1.0);
        assert!(parse_hap_names(&self.hap_names).is_some());
    }
}

//TODO use PathBuf
#[derive(clap::Args, Debug)]
pub struct TrioSettings {
    /// GFA file
    #[clap(short, long)]
    graph: PathBuf,

    /// Parental markers file
    #[clap(short, long)]
    markers: PathBuf,

    /// Marker-based annotation output file
    #[clap(long)]
    init_assign: Option<PathBuf>,

    /// Refined annotation output file
    #[clap(long)]
    refined_assign: Option<PathBuf>,

    /// Final annotation output file
    #[clap(long)]
    final_assign: Option<PathBuf>,

    /// Comma separated haplotype names to be used in outputs (default: "mat,pat")
    #[clap(long, default_value_t = String::from("mat,pat"))]
    hap_names: String,

    /// Marker-assisted extracted haplo-paths
    #[clap(long, short)]
    paths: Option<PathBuf>,

    /// Use GAF ([<>]<name1>)+ format for paths
    #[clap(long)]
    gaf_format: bool,

    /// Minimal number of parent-specific markers required for assigning parental group to a node
    #[clap(long, default_value_t = 10)]
    marker_cnt: usize,

    /// Require at least (node_length / <value>) markers within the node for parental group assignment
    #[clap(long, default_value_t = 10_000)]
    marker_sparsity: usize,

    /// Sets minimal marker excess for assigning a parental group to <value>:1
    #[clap(long, default_value_t = 5.0)]
    marker_ratio: f64,

    /// Longer nodes are unlikely to be spurious and likely to be reliably assigned based on markers (used in HOMOZYGOUS node labeling)
    #[clap(long, default_value_t = 200_000)]
    trusted_len: usize,

    /// Nodes with coverage below <coeff> * <weighted mean coverage of 'solid' nodes> can not be 'reclassified' as homozygous.
    /// Negative turns off reclassification, 0. disables coverage check
    #[clap(long, default_value_t = 1.5)]
    suspect_homozygous_cov_coeff: f64,

    /// Longer nodes can not be classified as homozygous
    #[clap(long, default_value_t = 2_000_000)]
    max_homozygous_len: usize,

    //TODO maybe check that it is > trusted_len
    /// Longer nodes are unlikely to represent repeats, polymorphic variants, etc (used to seed and guide the path search)
    #[clap(long, default_value_t = 500_000)]
    solid_len: usize,

    /// Sets minimal marker excess for assigning a parental group of solid nodes to <value>:1.
    /// Must be <= marker_ratio (by default == marker_ratio)
    #[clap(long)]
    solid_ratio: Option<f64>,

    /// Solid nodes with coverage below <coeff> * <weighted mean coverage of 'solid' nodes> can not be classified as homozygous.
    /// 0. disables check
    #[clap(long, default_value_t = 1.5)]
    solid_homozygous_cov_coeff: f64,

    /// Minimal node length for assigning ISSUE label
    #[clap(long, default_value_t = 50_000)]
    issue_len: usize,

    /// Minimal number of markers for assigning ISSUE label (by default == marker_cnt, will typically be set to a value >= marker_cnt)
    #[clap(long)]
    issue_cnt: Option<usize>,

    /// Require at least (node_length / <value>) markers for assigning ISSUE label (by default == marker_sparsity, will typically be set to a value >= marker_sparsity)
    #[clap(long)]
    issue_sparsity: Option<usize>,

    /// Require primary marker excess BELOW <value>:1 for assigning ISSUE label. Must be <= marker_ratio (by default == marker_ratio)
    #[clap(long)]
    issue_ratio: Option<f64>,

    /// Try to fill in small ambiguous bubbles
    #[clap(long)]
    try_fill_bubbles: bool,

    /// Do not fill bubble if source or sink is non-solid, non-homozygous and has coverage above <coeff> * <weighted mean coverage of 'solid' nodes>.
    /// Negative disables check, 0. makes it fail
    #[clap(long, default_value_t = 1.5)]
    max_unique_cov_coeff: f64,

    /// Bubbles including a longer alternative sequence will not be filled
    #[clap(long, default_value_t = 50_000)]
    fillable_bubble_len: usize,

    /// Bubbles with bigger difference between alternatives' lengths will not be filled
    #[clap(long, default_value_t = 200)]
    fillable_bubble_diff: usize,

    /// Heterozygous bubbles including a longer alternative sequence will not be filled (by default equal to fillable_bubble_len)
    #[clap(long)]
    het_fill_bubble_len: Option<usize>,

    /// Heterozygous bubbles with bigger difference between alternatives' lengths will not be filled (by default equal to fillable_bubble_diff)
    #[clap(long)]
    het_fill_bubble_diff: Option<usize>,

    /// During bubble filling ignore simple sides of bubbles with coverage less than source/sink average divided by this value
    /// 0. disables check
    #[clap(long, default_value_t = 5.0)]
    good_side_cov_gap: f64,

    /// Minimal introducible gap size (number of Ns reported). If the gap size estimate is smaller it will be artificially increased to this value.
    #[clap(long, default_value_t = 1000)]
    min_gap_size: usize,

    /// Default gap size, which will be output in cases where reasonable estimate is not possible or (more likely) hasn't been implemented yet.
    #[clap(long, default_value_t = 5000)]
    default_gap_size: usize,

    /// Assign tangles flanked by solid nodes from the same class
    #[clap(long)]
    assign_tangles: bool,

    /// Allow dead-end nodes in the tangles
    #[clap(long)]
    tangle_allow_deadend: bool,

    /// Check that inner tangle nodes are either unassigned or assigned to correct class
    #[clap(long)]
    tangle_check_inner: bool,

    /// Prevent reassignment of nodes
    #[clap(long)]
    tangle_prevent_reassign: bool,
}

impl TrioSettings {
    pub fn validate(&self) {
        if let Some(issue_ratio) = self.issue_ratio {
            assert!(
                issue_ratio <= self.marker_ratio,
                "--issue-ratio can't be set to a value higher than --marker-ratio"
            );
        }

        if let Some(solid_ratio) = self.solid_ratio {
            assert!(
                solid_ratio <= self.marker_ratio,
                "--solid-ratio can't be set to a value higher than --marker-ratio"
            );

            if solid_ratio < self.issue_ratio.unwrap_or(self.marker_ratio) {
                warn!(
                    "Specified --solid-ratio value is smaller than --issue-ratio. \
                    Please double-check the logic and consider specifying smaller --issue-ratio."
                );
            }
        }

        assert!(self.good_side_cov_gap >= 0.);
        assert!(self.solid_homozygous_cov_coeff >= 0.);
    }
}

fn read_graph(graph_fn: &PathBuf) -> Result<Graph, Box<dyn Error>> {
    info!("Reading graph from {}", graph_fn.to_str().unwrap());
    let g = Graph::read_sanitize(&fs::read_to_string(graph_fn)?);

    info!("Graph read successfully");
    info!("Node count: {}", g.node_cnt());
    info!("Link count: {}", g.link_cnt());
    Ok(g)
}

fn output_coloring(
    g: &Graph,
    assignments: &trio::AssignmentStorage,
    file_name: &PathBuf,
    hap_names: &(&str, &str),
) -> Result<(), std::io::Error> {
    let mut output = BufWriter::new(File::create(file_name)?);
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
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}",
                n.name,
                group_str(Some(assign.group), hap_names).to_uppercase(),
                n.length,
                assign.info,
                color
            )?;
        }
    }
    Ok(())
}

pub fn augment_by_path_search(
    g: &Graph,
    assignments: trio::AssignmentStorage,
    settings: HaploSearchSettings,
) -> trio::AssignmentStorage {
    info!("Augmenting node annotation by path search. Round 1.");
    let assignments = augment_by_path_search_round(g, assignments, settings);
    info!("Augmenting node annotation by path search. Round 2.");
    augment_by_path_search_round(g, assignments, settings)
}

fn augment_by_path_search_round(
    g: &Graph,
    assignments: trio::AssignmentStorage,
    settings: HaploSearchSettings,
) -> trio::AssignmentStorage {
    let mut path_searcher =
        HaploSearcher::new(g, &assignments, settings.assigning_stage_adjusted(), None);

    path_searcher.find_all();
    let node_usage = path_searcher.take_used();
    augment_assignments(g, assignments, &node_usage, true)
}

fn augment_assignments(
    g: &Graph,
    mut assignments: trio::AssignmentStorage,
    extra_assignments: &trio::AssignmentStorage,
    exclude_homozygous: bool,
) -> trio::AssignmentStorage {
    for node_id in extra_assignments.assigned() {
        let tentative_group = extra_assignments.group(node_id).unwrap();
        assert!(tentative_group != TrioGroup::ISSUE);
        //any mixed assignment has chance to be erroneous due to graph issues
        if exclude_homozygous && !tentative_group.is_definite() {
            continue;
        }
        match assignments.group(node_id) {
            None => {
                debug!(
                    "Assigning tentative group {:?} to node {}",
                    tentative_group,
                    g.name(node_id)
                );
                assignments.assign(node_id, tentative_group, "PathSearch");
            }
            Some(init_group) => {
                assert!(init_group == tentative_group || init_group == trio::TrioGroup::HOMOZYGOUS)
            }
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

fn parse_hap_names(hap_names_s: &str) -> Option<(&str, &str)> {
    let mut split = hap_names_s.split(',');
    Some((split.next()?, split.next()?))
}

fn group_str<'a>(o_g: Option<TrioGroup>, hap_names: &'a (&'a str, &'a str)) -> &'a str {
    match o_g {
        Some(TrioGroup::MATERNAL) => hap_names.0,
        Some(TrioGroup::PATERNAL) => hap_names.1,
        Some(TrioGroup::HOMOZYGOUS) => "hom",
        Some(TrioGroup::ISSUE) => "issue",
        _ => "na",
    }
}

pub fn write_paths(
    g: &Graph,
    haplo_paths: Vec<trio_walk::HaploPath>,
    assignments: &trio::AssignmentStorage,
    node_usage: &trio::AssignmentStorage,
    output: &PathBuf,
    gaf_format: bool,
    hap_names: &(&str, &str),
) -> Result<(), std::io::Error> {
    //FIXME buffer
    let mut output = File::create(output)?;
    writeln!(output, "name\tpath\tassignment")?;
    for (path, node_id, group) in haplo_paths {
        assert!(path.vertices().contains(&Vertex::forward(node_id)));
        //info!("Identified {:?} path: {}", group, path.print(&g));
        writeln!(
            output,
            "{}_from_{}\t{}\t{}",
            group_str(Some(group), hap_names),
            g.node(node_id).name,
            path.print_format(g, gaf_format),
            group_str(Some(group), hap_names).to_uppercase()
        )?;
    }

    let mut write_node = |n: &Node, group: Option<TrioGroup>| {
        writeln!(
            output,
            "{}_unused_{}\t{}\t{}",
            group_str(group, hap_names),
            n.name,
            Direction::format_node(&n.name, Direction::FORWARD, gaf_format),
            group_str(group, hap_names).to_uppercase()
        )
    };

    for (node_id, n) in g.all_nodes().enumerate() {
        let haplopath_assign = node_usage.group(node_id);
        match assignments.group(node_id) {
            None | Some(TrioGroup::ISSUE) => {
                assert!(!node_usage.contains(node_id));
                debug!(
                    "Node: {} length: {} not assigned to any haplotype (adding trivial NA path)",
                    n.name, n.length
                );
                write_node(g.node(node_id), None)?;
            }
            Some(assign) => {
                if TrioGroup::compatible(assign, TrioGroup::MATERNAL)
                    //not present in haplopaths paths or incompatible
                    && haplopath_assign.map_or(true,
                        |x| TrioGroup::incompatible(x, TrioGroup::MATERNAL))
                {
                    debug!("Node: {} length: {} not present in MATERNAL haplo-paths (adding trivial MATERNAL path)",
                        n.name, n.length);
                    write_node(g.node(node_id), Some(TrioGroup::MATERNAL))?;
                }
                if TrioGroup::compatible(assign, TrioGroup::PATERNAL)
                    //not present in haplopaths paths or incompatible
                    && haplopath_assign.map_or(true,
                        |x| TrioGroup::incompatible(x, TrioGroup::PATERNAL))
                {
                    debug!("Node: {} length: {} not present in PATERNAL haplo-paths (adding trivial PATERNAL path)",
                        n.name, n.length);
                    write_node(g.node(node_id), Some(TrioGroup::PATERNAL))?;
                }
            }
        }
    }
    Ok(())
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

    let hap_names =
        parse_hap_names(&settings.hap_names).expect("Problem while parsing haplotype names");

    info!(
        "Reading trio marker information from {}",
        &settings.markers.to_str().unwrap()
    );
    let trio_infos = trio::read_trio(&settings.markers)?;

    let solid_cov_est = weighted_mean_solid_cov(&g, settings.solid_len);
    if settings.suspect_homozygous_cov_coeff > 0. || settings.solid_homozygous_cov_coeff > 0. {
        info!("Coverage estimate based on long nodes was {solid_cov_est}");
        if solid_cov_est == 0. {
            warn!("Looks like the graph didn't have coverage information, which we were hoping to use. \
                    Consider providing it or changing --suspect-homozygous-cov-coeff and --solid-homozygous-cov-coeff");
        }
    }

    let suspect_homozygous_cov = if settings.suspect_homozygous_cov_coeff < 0. {
        None
    } else {
        Some(settings.suspect_homozygous_cov_coeff * solid_cov_est)
    };

    let solid_homozygous_cov = settings.solid_homozygous_cov_coeff * solid_cov_est;

    info!("Assigning initial parental groups to the nodes");
    let assignments = trio::assign_parental_groups(
        &g,
        &trio_infos,
        &GroupAssignmentSettings {
            assign_cnt: settings.marker_cnt,
            assign_sparsity: settings.marker_sparsity,
            assign_ratio: settings.marker_ratio,
            solid_ratio: settings.solid_ratio.unwrap_or(settings.marker_ratio),
            issue_len: settings.issue_len,
            issue_cnt: settings.issue_cnt.unwrap_or(settings.marker_cnt),
            issue_sparsity: settings.issue_sparsity.unwrap_or(settings.marker_sparsity),
            issue_ratio: settings.issue_ratio.unwrap_or(settings.marker_ratio),
        },
        settings.solid_len,
        solid_homozygous_cov,
    );

    let raw_cnts = trio_infos
        .into_iter()
        .map(|ti| (g.name2id(&ti.node_name), ti))
        .collect::<HashMap<usize, trio::TrioInfo>>();

    if let Some(output) = &settings.init_assign {
        info!(
            "Writing initial node annotation to {}",
            output.to_str().unwrap()
        );
        output_coloring(&g, &assignments, output, &hap_names)?;
    }

    info!("Marking homozygous nodes");
    let assigner = trio::HomozygousAssigner::new(
        &g,
        assignments,
        settings.trusted_len,
        suspect_homozygous_cov,
        settings.solid_len,
        solid_homozygous_cov,
        settings.max_homozygous_len,
    );

    let assignments = assigner.run();

    let mut search_settings = HaploSearchSettings {
        solid_len: settings.solid_len,
        trusted_len: settings.trusted_len,
        fill_bubbles: settings.try_fill_bubbles,
        fillable_bubble_len: settings.fillable_bubble_len,
        fillable_bubble_diff: settings.fillable_bubble_diff,
        het_fill_bubble_len: settings
            .het_fill_bubble_len
            .unwrap_or(settings.fillable_bubble_len),
        het_fill_bubble_diff: settings
            .het_fill_bubble_diff
            .unwrap_or(settings.fillable_bubble_diff),
        good_side_cov_gap: settings.good_side_cov_gap,
        min_gap_size: settings.min_gap_size as i64,
        default_gap_size: settings.default_gap_size as i64,
        ..HaploSearchSettings::default()
    };

    if search_settings.fill_bubbles {
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
        info!(
            "Maximal 'unique' coverage for bubble filling set to {}",
            search_settings.max_unique_cov
        );
        if search_settings.max_unique_cov == 0. {
            info!("Will only fill bubbles between solid or homozygous nodes");
        }
    }

    let assignments = augment_by_path_search(&g, assignments, search_settings);

    let assignments = if settings.assign_tangles {
        assign_short_node_tangles(
            &g,
            assignments,
            settings.solid_len,
            TangleAssignmentSettings {
                allow_deadend: settings.tangle_allow_deadend,
                check_inner: settings.tangle_check_inner,
                allow_reassign: !settings.tangle_prevent_reassign,
            },
        )
    } else {
        assignments
    };

    if let Some(output) = &settings.refined_assign {
        info!(
            "Writing refined node annotation to {}",
            output.to_str().unwrap()
        );
        output_coloring(&g, &assignments, output, &hap_names)?;
    }
    let mut path_searcher = HaploSearcher::new(&g, &assignments, search_settings, Some(&raw_cnts));

    let haplo_paths = path_searcher.find_all();
    let node_usage = path_searcher.take_used();

    let assignments = augment_assignments(&g, assignments, &node_usage, false);

    if let Some(output) = &settings.final_assign {
        info!(
            "Writing final node annotation to {}",
            output.to_str().unwrap()
        );
        output_coloring(&g, &assignments, output, &hap_names)?;
    }

    if let Some(output) = &settings.paths {
        info!("Outputting haplo-paths to {}", output.to_str().unwrap());
        write_paths(
            &g,
            haplo_paths,
            &assignments,
            &node_usage,
            output,
            settings.gaf_format,
            &hap_names,
        )?;
    }

    info!("All done");
    Ok(())
}

/// Length-weighted coverage is robust to graph fragmentation: splitting a
/// sequence into many nodes does not give that sequence extra influence.
fn weighted_median_cov(g: &Graph, node_ids: impl Iterator<Item = usize>) -> Option<f64> {
    let mut cov_len: Vec<(f64, usize)> = node_ids
        .map(|node_id| (g.node(node_id).coverage, g.node(node_id).length))
        .filter(|(coverage, length)| *coverage > 0.0 && *length > 0)
        .collect();
    if cov_len.is_empty() {
        return None;
    }
    cov_len.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let total_len: usize = cov_len.iter().map(|(_, length)| length).sum();
    let mut seen = 0;
    for (coverage, length) in cov_len {
        seen += length;
        if seen * 2 >= total_len {
            return Some(coverage);
        }
    }
    unreachable!()
}

fn dual_block_hash(g: &Graph, block: &pseudo_hap::LinearBlock, seed: u64) -> u64 {
    // FNV-1a over sorted node names makes the decision independent of graph
    // traversal direction and stable across platforms.
    let mut names: Vec<&str> = block.all_nodes().map(|node_id| g.name(node_id)).collect();
    names.sort_unstable();
    let mut hash = 0xcbf29ce484222325u64 ^ seed;
    for name in names {
        for byte in name.as_bytes().iter().chain(std::iter::once(&0u8)) {
            hash ^= *byte as u64;
            hash = hash.wrapping_mul(0x100000001b3);
        }
    }
    hash
}

/// Find connected graph components not covered by a pseudo-hap bubble block.
/// Sex chromosomes and other haploid components commonly enter this path.
fn unassigned_components(g: &Graph, assignments: &trio::AssignmentStorage) -> Vec<Vec<usize>> {
    let mut neighbors = vec![Vec::<usize>::new(); g.node_cnt()];
    for link in g.all_links() {
        if link.start.node_id != link.end.node_id {
            neighbors[link.start.node_id].push(link.end.node_id);
            neighbors[link.end.node_id].push(link.start.node_id);
        }
    }

    let mut seen = vec![false; g.node_cnt()];
    let mut components = Vec::new();
    for start in 0..g.node_cnt() {
        if seen[start] || assignments.contains(start) {
            continue;
        }
        let mut stack = vec![start];
        let mut component = Vec::new();
        seen[start] = true;
        while let Some(node_id) = stack.pop() {
            component.push(node_id);
            for &next in &neighbors[node_id] {
                if !seen[next] && !assignments.contains(next) {
                    seen[next] = true;
                    stack.push(next);
                }
            }
        }
        components.push(component);
    }
    components
}

fn undirected_component_ids(g: &Graph) -> Vec<usize> {
    let mut neighbors = vec![Vec::<usize>::new(); g.node_cnt()];
    for link in g.all_links() {
        if link.start.node_id != link.end.node_id {
            neighbors[link.start.node_id].push(link.end.node_id);
            neighbors[link.end.node_id].push(link.start.node_id);
        }
    }
    let mut component_of = vec![usize::MAX; g.node_cnt()];
    let mut component_id = 0;
    for start in 0..g.node_cnt() {
        if component_of[start] != usize::MAX {
            continue;
        }
        component_of[start] = component_id;
        let mut stack = vec![start];
        while let Some(node_id) = stack.pop() {
            for &next in &neighbors[node_id] {
                if component_of[next] == usize::MAX {
                    component_of[next] = component_id;
                    stack.push(next);
                }
            }
        }
        component_id += 1;
    }
    component_of
}

/// Return the distinct unresolved neighbors reachable from one oriented side
/// of a node.  GFA links are bidirected, so the forward and reverse vertices
/// represent the two physical sides of the sequence.
fn unresolved_side_neighbors(
    g: &Graph,
    assignments: &trio::AssignmentStorage,
    vertex: Vertex,
) -> HashSet<usize> {
    g.outgoing_edges(vertex)
        .into_iter()
        .map(|link| link.end.node_id)
        .filter(|&node_id| node_id != vertex.node_id && !assignments.contains(node_id))
        .collect()
}

fn side_neighbors(g: &Graph, vertex: Vertex) -> HashSet<usize> {
    g.outgoing_edges(vertex)
        .into_iter()
        .map(|link| link.end.node_id)
        .filter(|&node_id| node_id != vertex.node_id)
        .collect()
}

/// Detect collapsed junctions missed by pseudo-hap decomposition.
///
/// Coverage alone is insufficient: a high-copy linear or repetitive node
/// should not split a component.  Require at least two unresolved branches on
/// each oriented side and haploid-depth support on both sides.  This captures
/// terminal "bowtie" structures while leaving recognized bubbles untouched.
fn unresolved_shared_anchors(
    g: &Graph,
    assignments: &trio::AssignmentStorage,
    haploid_cov: f64,
    shared_cov_coeff: f64,
) -> Vec<usize> {
    if haploid_cov <= 0.0 {
        return Vec::new();
    }
    let max_single_cov = shared_cov_coeff * haploid_cov;
    let mut anchors = Vec::new();
    for node_id in 0..g.node_cnt() {
        if assignments.contains(node_id) || g.node(node_id).coverage < max_single_cov {
            continue;
        }
        let forward = unresolved_side_neighbors(g, assignments, Vertex::forward(node_id));
        let reverse = unresolved_side_neighbors(g, assignments, Vertex::reverse(node_id));
        // Resolve only binary choices.  Higher-order repetitive junctions do
        // not have a unique complementary split without additional evidence.
        if forward.len() != 2 || reverse.len() != 2 {
            continue;
        }
        let all_single_copy_neighbors = |neighbors: &HashSet<usize>| {
            neighbors.iter().all(|&neighbor| {
                let coverage = g.node(neighbor).coverage;
                coverage > 0.0 && coverage < max_single_cov
            })
        };
        // All four immediate arms must look single-copy.  Allowing only one
        // supporting arm per side admits high-copy/zero-coverage rDNA tangles
        // whose removal can repartition several otherwise-correct components.
        if all_single_copy_neighbors(&forward) && all_single_copy_neighbors(&reverse) {
            anchors.push(node_id);
        }
    }
    anchors
}

fn stable_node_hash(g: &Graph, node_ids: &[usize], seed: u64) -> u64 {
    let mut names: Vec<&str> = node_ids.iter().map(|&node_id| g.name(node_id)).collect();
    names.sort_unstable();
    let mut hash = 0xcbf29ce484222325u64 ^ seed;
    for name in names {
        for byte in name.as_bytes().iter().chain(std::iter::once(&0u8)) {
            hash ^= *byte as u64;
            hash = hash.wrapping_mul(0x100000001b3);
        }
    }
    hash
}

fn read_dual_homology(
    g: &Graph,
    path: &PathBuf,
) -> Result<Vec<(usize, usize, usize)>, Box<dyn Error>> {
    let name_to_id: HashMap<&str, usize> = (0..g.node_cnt())
        .map(|node_id| (g.name(node_id), node_id))
        .collect();
    let mut pairs = Vec::new();
    for (line_number, line) in fs::read_to_string(path)?.lines().enumerate() {
        if line_number == 0 && line.starts_with("node1\t") {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.is_empty() {
            continue;
        }
        if fields.len() != 3 {
            return Err(format!(
                "{}:{}: expected node1, node2, and homology length",
                path.display(),
                line_number + 1
            )
            .into());
        }
        let Some(&first) = name_to_id.get(fields[0]) else {
            warn!("Ignoring homology for absent graph node {}", fields[0]);
            continue;
        };
        let Some(&second) = name_to_id.get(fields[1]) else {
            warn!("Ignoring homology for absent graph node {}", fields[1]);
            continue;
        };
        let homology_length = fields[2].parse::<usize>().map_err(|error| {
            format!(
                "{}:{}: invalid homology length: {error}",
                path.display(),
                line_number + 1
            )
        })?;
        if first != second && homology_length > 0 {
            pairs.push((first, second, homology_length));
        }
    }
    Ok(pairs)
}

/// Assign the two single-copy branches on each side of a validated collapsed
/// anchor to opposite outputs.  These are topology constraints, not phasing
/// evidence: existing colors select the orientation that keeps the most
/// sequence unchanged, with the seed used only for exact ties.
fn assign_binary_choices(
    g: &Graph,
    anchors: &[usize],
    assignments: &mut trio::AssignmentStorage,
    seed: u64,
    hap_lengths: &mut [usize; 2],
    max_single_cov: f64,
) {
    let mut constraints = vec![Vec::<(usize, bool)>::new(); g.node_cnt()];
    for &node_id in anchors {
        for vertex in [Vertex::forward(node_id), Vertex::reverse(node_id)] {
            let choices = side_neighbors(g, vertex);
            if choices.len() == 2 {
                let mut pair: Vec<usize> = choices.into_iter().collect();
                pair.sort_by_key(|&node_id| g.name(node_id));
                let first = pair[0];
                let second = pair[1];
                let single_copy = |choice: usize| {
                    let coverage = g.node(choice).coverage;
                    coverage > 0.0 && coverage < max_single_cov
                };
                if first != second && single_copy(first) && single_copy(second) {
                    constraints[first].push((second, true));
                    constraints[second].push((first, true));
                }
            }
        }
    }

    let mut colors = vec![None::<bool>; g.node_cnt()];
    for start in 0..g.node_cnt() {
        if constraints[start].is_empty() || colors[start].is_some() {
            continue;
        }
        let mut queue = VecDeque::from([start]);
        let mut members = Vec::new();
        let mut conflict = false;
        colors[start] = Some(false);
        while let Some(component_id) = queue.pop_front() {
            members.push(component_id);
            let color = colors[component_id].unwrap();
            for &(next, opposite) in &constraints[component_id] {
                let expected = color ^ opposite;
                match colors[next] {
                    None => {
                        colors[next] = Some(expected);
                        queue.push_back(next);
                    }
                    Some(next_color) if next_color != expected => conflict = true,
                    Some(_) => {}
                }
            }
        }
        if conflict {
            warn!(
                "Contradictory dual-assembly choices across {} residual segments; leaving them for the conservative fallback",
                members.len()
            );
            for &component_id in &members {
                colors[component_id] = None;
            }
            continue;
        }

        let agreement_score = |flip: bool| {
            members
                .iter()
                .map(|&node_id| {
                    let desired = colors[node_id].unwrap() ^ flip;
                    let agrees = matches!(
                        (assignments.group(node_id), desired),
                        (Some(TrioGroup::MATERNAL), false) | (Some(TrioGroup::PATERNAL), true)
                    );
                    if agrees {
                        g.node_length(node_id)
                    } else {
                        0
                    }
                })
                .sum::<usize>()
        };
        let keep_score = agreement_score(false);
        let flip_score = agreement_score(true);
        let flip = if keep_score == flip_score {
            stable_node_hash(g, &members, seed) & 1 == 1
        } else {
            flip_score > keep_score
        };
        for &node_id in &members {
            let group_idx = (colors[node_id].unwrap() ^ flip) as usize;
            let group = if group_idx == 0 {
                TrioGroup::MATERNAL
            } else {
                TrioGroup::PATERNAL
            };
            assignments.assign(node_id, group, "dual-binary-choice");
            hap_lengths[group_idx] += g.node_length(node_id);
        }
        info!(
            "Assigned {} nodes using complementary binary graph choices",
            members.len()
        );
    }
}

/// Reproduce the original component-level fallback without mutating the real
/// assignments.  New shared-anchor choices can then replace only the affected
/// components while every unrelated fallback decision remains stable.
fn baseline_component_groups(
    g: &Graph,
    assignments: &trio::AssignmentStorage,
    haploid_cov: f64,
    haploid_component_cov_coeff: f64,
    initial_hap_lengths: [usize; 2],
    seed: u64,
) -> Vec<Option<TrioGroup>> {
    let components = unassigned_components(g, assignments);
    let mut component_order: Vec<usize> = (0..components.len()).collect();
    component_order.sort_by_key(|&component_id| {
        usize::MAX
            - components[component_id]
                .iter()
                .map(|&node_id| g.node_length(node_id))
                .sum::<usize>()
    });

    let mut groups = vec![None; g.node_cnt()];
    let mut hap_lengths = initial_hap_lengths;
    for component_id in component_order {
        let component = &components[component_id];
        let component_len = component
            .iter()
            .map(|&node_id| g.node_length(node_id))
            .sum::<usize>();
        let component_cov = weighted_median_cov(g, component.iter().copied()).unwrap_or(0.0);
        let group = if haploid_cov > 0.0
            && component_cov > 0.0
            && component_cov < haploid_component_cov_coeff * haploid_cov
        {
            let group_idx = if hap_lengths[0] == hap_lengths[1] {
                let mut hash = seed ^ 0x9e3779b97f4a7c15;
                let mut names: Vec<&str> = component.iter().map(|&id| g.name(id)).collect();
                names.sort_unstable();
                for name in names {
                    for byte in name.as_bytes() {
                        hash = hash.rotate_left(5) ^ (*byte as u64);
                    }
                }
                (hash & 1) as usize
            } else if hap_lengths[0] < hap_lengths[1] {
                0
            } else {
                1
            };
            hap_lengths[group_idx] += component_len;
            Some(if group_idx == 0 {
                TrioGroup::MATERNAL
            } else {
                TrioGroup::PATERNAL
            })
        } else if haploid_cov > 0.0
            && component_cov >= haploid_component_cov_coeff * haploid_cov
            && component_cov <= 2.5 * haploid_cov
        {
            Some(TrioGroup::HOMOZYGOUS)
        } else {
            None
        };
        for &node_id in component {
            groups[node_id] = group;
        }
    }
    groups
}

/// Build two complementary assemblies without external phasing evidence.
///
/// Bubble interiors are split between outputs, shared boundaries are marked
/// homozygous, and untouched components are classified by their coverage.
pub fn run_dual_analysis(settings: &DualSettings) -> Result<(), Box<dyn Error>> {
    let g = read_graph(&settings.graph)?;
    let hap_names =
        parse_hap_names(&settings.hap_names).expect("Problem while parsing haplotype names");
    let blocks = pseudo_hap::pseudo_hap_decompose(&g, settings.unique_block_len);
    for block in &blocks {
        debug!(
            "Dual pseudo-hap block primary=[{}] alternate=[{}] shared=[{}]",
            block
                .primary_exclusive_nodes()
                .map(|node_id| g.name(node_id))
                .collect::<Vec<_>>()
                .join(","),
            block
                .known_alt_nodes()
                .iter()
                .map(|&node_id| g.name(node_id))
                .collect::<Vec<_>>()
                .join(","),
            block
                .shared_nodes()
                .iter()
                .map(|&node_id| g.name(node_id))
                .collect::<Vec<_>>()
                .join(",")
        );
    }

    // Bubble sides represent single-copy alternatives even when most of the
    // graph is homozygous, so they are the preferred haploid-depth estimator.
    let bubble_nodes: HashSet<usize> = blocks
        .iter()
        .filter(|block| !block.known_alt_nodes().is_empty())
        .flat_map(|block| {
            block
                .primary_exclusive_nodes()
                .chain(block.known_alt_nodes().iter().copied())
        })
        .collect();
    let haploid_cov = if let Some(coverage) = weighted_median_cov(&g, bubble_nodes.iter().copied())
    {
        info!("Estimated haploid coverage from pseudo-hap bubbles: {coverage:.3}");
        coverage
    } else {
        let graph_median = weighted_median_cov(&g, 0..g.node_cnt()).unwrap_or(0.0);
        if settings.haploid {
            warn!("No covered pseudo-hap bubbles were available for estimating haploid coverage; --haploid is set, so using the graph-wide weighted median.");
            graph_median
        } else {
            warn!("No covered pseudo-hap bubbles were available for estimating haploid coverage; using half the graph-wide weighted median.");
            graph_median / 2.0
        }
    };
    info!("Dual-mode haploid coverage is {haploid_cov:.3}");

    let mut assignments = trio::AssignmentStorage::new();
    let mut hap_lengths = [0usize, 0usize];
    let component_of = undirected_component_ids(&g);
    let mut linear_blocks = Vec::<(Vec<usize>, usize)>::new();

    for block in &blocks {
        if block.known_alt_nodes().is_empty() {
            // A long linear block is still informative in dual mode even
            // without an explicit ALT path.  Single-copy blocks represent a
            // haploid component/arm and belong in one output; diploid-depth
            // blocks are collapsed sequence shared by both outputs.  Skipping
            // these blocks lets a later connected-component fallback flatten
            // several distinct chromosome arms to one color.
            let block_nodes: Vec<usize> = block.primary_exclusive_nodes().collect();
            let block_cov = weighted_median_cov(&g, block_nodes.iter().copied()).unwrap_or(0.0);
            let max_single_cov = settings.haploid_component_cov_coeff * haploid_cov;
            if block_cov > 0.0 && block_cov < max_single_cov {
                let block_len = block_nodes
                    .iter()
                    .map(|&node_id| g.node_length(node_id))
                    .sum();
                linear_blocks.push((block_nodes, block_len));
            } else if block_cov >= max_single_cov && block_cov <= 2.5 * haploid_cov {
                for node_id in block_nodes {
                    assignments.assign(node_id, TrioGroup::HOMOZYGOUS, "dual-linear-shared");
                }
            }
            continue;
        }
        // Swap primary/ALT independently per block.  Hashing sorted node names
        // makes the split reproducible without favoring traversal order.
        let flip = dual_block_hash(&g, block, settings.seed) & 1 == 1;
        let primary_group = if flip {
            TrioGroup::PATERNAL
        } else {
            TrioGroup::MATERNAL
        };
        let alt_group = if flip {
            TrioGroup::MATERNAL
        } else {
            TrioGroup::PATERNAL
        };

        for node_id in block.shared_nodes() {
            assignments.update_group(*node_id, TrioGroup::HOMOZYGOUS);
        }
        for node_id in block.primary_exclusive_nodes() {
            assignments.update_group(node_id, primary_group);
            hap_lengths[if primary_group == TrioGroup::MATERNAL {
                0
            } else {
                1
            }] += g.node_length(node_id);
        }
        for &node_id in block.known_alt_nodes() {
            assignments.update_group(node_id, alt_group);
            hap_lengths[if alt_group == TrioGroup::MATERNAL {
                0
            } else {
                1
            }] += g.node_length(node_id);
        }
    }

    let mut paired_with = vec![None::<usize>; linear_blocks.len()];
    let mut mashmap_pair = HashSet::<(usize, usize)>::new();
    if let Some(homology_path) = &settings.homology {
        // Convert strict node homologies into block-level scores.  The Python
        // selector has already applied the Hi-C pipeline's interval coverage,
        // adjacency, and clear-best filters; aggregation here handles the
        // uncommon case where a pseudo-hap block contains several nodes.
        let mut linear_of_node = vec![usize::MAX; g.node_cnt()];
        for (linear_id, (nodes, _)) in linear_blocks.iter().enumerate() {
            for &node_id in nodes {
                linear_of_node[node_id] = linear_id;
            }
        }
        let mut scores = HashMap::<(usize, usize), usize>::new();
        for (first_node, second_node, homology_length) in read_dual_homology(&g, homology_path)? {
            let first = linear_of_node[first_node];
            let second = linear_of_node[second_node];
            if first == usize::MAX || second == usize::MAX || first == second {
                continue;
            }
            let pair = if first < second {
                (first, second)
            } else {
                (second, first)
            };
            *scores.entry(pair).or_default() += homology_length;
        }

        // Reapply the established 1.5x clear-best rule at block level, then
        // greedily take the strongest non-conflicting pairs.  This prevents a
        // repeat-rich block from being paired with several alternatives.
        let mut best_two = vec![[0usize; 2]; linear_blocks.len()];
        for (&(first, second), &score) in &scores {
            for linear_id in [first, second] {
                if score > best_two[linear_id][0] {
                    best_two[linear_id][1] = best_two[linear_id][0];
                    best_two[linear_id][0] = score;
                } else if score > best_two[linear_id][1] {
                    best_two[linear_id][1] = score;
                }
            }
        }
        let mut candidates: Vec<((usize, usize), usize)> = scores.into_iter().collect();
        candidates.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
        for ((first, second), score) in candidates {
            let clear_for = |linear_id: usize| {
                score == best_two[linear_id][0]
                    && (best_two[linear_id][1] == 0
                        || score.saturating_mul(2) >= best_two[linear_id][1].saturating_mul(3))
            };
            if paired_with[first].is_none()
                && paired_with[second].is_none()
                && (clear_for(first) || clear_for(second))
            {
                paired_with[first] = Some(second);
                paired_with[second] = Some(first);
                mashmap_pair.insert((first, second));
                info!(
                    "Paired dual linear blocks with {} bp of strict MashMap homology",
                    score
                );
            }
        }
    } else {
        warn!(
            "No --homology input was supplied; falling back to similar-length linear-block pairing"
        );
        // Retain a deterministic standalone fallback for direct Rukki use and
        // old saved graphs that no longer have sequence.  The Verkko workflow
        // always supplies MashMap evidence and does not enter this branch.
        let mut linear_by_component = HashMap::<usize, Vec<usize>>::new();
        for (linear_id, (nodes, _)) in linear_blocks.iter().enumerate() {
            if let Some(&node_id) = nodes.first() {
                linear_by_component
                    .entry(component_of[node_id])
                    .or_default()
                    .push(linear_id);
            }
        }
        for linear_ids in linear_by_component.values_mut() {
            linear_ids.sort_by_key(|&linear_id| usize::MAX - linear_blocks[linear_id].1);
            for pair in linear_ids.chunks_exact(2) {
                let first = pair[0];
                let second = pair[1];
                if linear_blocks[second].1.saturating_mul(5)
                    >= linear_blocks[first].1.saturating_mul(4)
                {
                    paired_with[first] = Some(second);
                    paired_with[second] = Some(first);
                }
            }
        }
    }

    let mut colored_linear = vec![false; linear_blocks.len()];
    for first in 0..linear_blocks.len() {
        if colored_linear[first] {
            continue;
        }
        let second = paired_with[first];
        let mut pair_nodes = linear_blocks[first].0.clone();
        if let Some(second) = second {
            pair_nodes.extend(linear_blocks[second].0.iter().copied());
        }
        let equal_length_flip = (stable_node_hash(&g, &pair_nodes, settings.seed) & 1) as usize;
        let first_group = if let Some(second) = second {
            // Haplotype 1 is paternal by convention.  For a divergent pair
            // such as X/Y, assigning the shorter member to hap1 places Y on
            // the paternal side; equal lengths retain seeded orientation.
            match linear_blocks[first].1.cmp(&linear_blocks[second].1) {
                std::cmp::Ordering::Less => 0,
                std::cmp::Ordering::Greater => 1,
                std::cmp::Ordering::Equal => equal_length_flip,
            }
        } else {
            equal_length_flip
        };
        let is_mashmap_pair = second
            .map(|second| mashmap_pair.contains(&(first.min(second), first.max(second))))
            .unwrap_or(false);
        for (linear_id, group_idx) in
            std::iter::once((first, first_group)).chain(second.map(|id| (id, first_group ^ 1)))
        {
            let group = if group_idx == 0 {
                TrioGroup::MATERNAL
            } else {
                TrioGroup::PATERNAL
            };
            let info = if is_mashmap_pair {
                "dual-homology-block"
            } else {
                "dual-linear-block"
            };
            for &node_id in &linear_blocks[linear_id].0 {
                assignments.assign(node_id, group, info);
                hap_lengths[group_idx] += g.node_length(node_id);
            }
            colored_linear[linear_id] = true;
        }
    }

    // Preserve the original fallback choices before introducing any new
    // anchors.  This prevents a repaired component from perturbing global
    // length balancing and flipping otherwise unrelated components.
    let baseline_groups = baseline_component_groups(
        &g,
        &assignments,
        haploid_cov,
        settings.haploid_component_cov_coeff,
        hap_lengths,
        settings.seed,
    );

    // Pseudo-hap decomposition intentionally leaves linear components alone.
    // Before the component-level fallback, identify high-coverage junctions
    // connecting haploid-depth branches.  Treating the entire mixed-copy
    // component as one unit would otherwise make all branches monochromatic.
    let shared_anchors = unresolved_shared_anchors(
        &g,
        &assignments,
        haploid_cov,
        settings.haploid_component_cov_coeff,
    );
    for &anchor in &shared_anchors {
        assignments.assign(anchor, TrioGroup::HOMOZYGOUS, "dual-shared-anchor");
        info!(
            "Collapsed shared anchor {} at {:.3}x",
            g.name(anchor),
            g.node(anchor).coverage
        );
    }

    // Assign binary branch choices around those anchors complementarily.
    // Remaining components retain the original copy-number fallback.
    // Restore the saved decision node by node.  Re-grouping after adding arm
    // seeds is unsafe because a new residual component can span several old
    // fallback groups through a repeat tangle; copying one group's color to
    // that whole component caused unrelated chromosomes to flip or collapse.
    for (node_id, baseline) in baseline_groups.iter().copied().enumerate() {
        if !assignments.contains(node_id) {
            if let Some(group) = baseline {
                assignments.assign(node_id, group, "dual-baseline");
            }
        }
    }

    assign_binary_choices(
        &g,
        &shared_anchors,
        &mut assignments,
        settings.seed,
        &mut hap_lengths,
        settings.haploid_component_cov_coeff * haploid_cov,
    );

    // Report only regions for which neither pseudo-hap decomposition, the
    // local anchor repair, nor the preserved copy-number fallback had a safe
    // assignment.
    for component in unassigned_components(&g, &assignments) {
        let component_len = component
            .iter()
            .map(|&node_id| g.node_length(node_id))
            .sum::<usize>();
        let component_cov = weighted_median_cov(&g, component.iter().copied()).unwrap_or(0.0);
        if component_len < settings.unique_block_len {
            debug!(
                "Short component of {component_len} bp at {component_cov:.3}x has uncertain copy number and remains unassigned"
            );
        } else {
            warn!(
                "Component of {component_len} bp at {component_cov:.3}x has uncertain copy number and remains unassigned"
            );
        }
    }

    let search_settings = HaploSearchSettings {
        max_unique_cov: settings.haploid_component_cov_coeff * haploid_cov,
        ..HaploSearchSettings::default()
    };
    let assignments = augment_by_path_search(&g, assignments, search_settings);
    let mut path_searcher = HaploSearcher::new(&g, &assignments, search_settings, None);
    let haplo_paths = path_searcher.find_all();
    let node_usage = path_searcher.take_used();
    let assignments = augment_assignments(&g, assignments, &node_usage, false);

    if let Some(output) = &settings.colors {
        output_coloring(&g, &assignments, output, &hap_names)?;
    }
    if let Some(output) = &settings.paths {
        write_paths(
            &g,
            haplo_paths,
            &assignments,
            &node_usage,
            output,
            settings.gaf_format,
            &hap_names,
        )?;
    }
    info!("All done");
    Ok(())
}

pub fn run_primary_alt_analysis(
    graph_fn: &PathBuf,
    colors_fn: &Option<String>,
    paths_fn: &Option<String>,
    gaf_paths: bool,
) -> Result<(), Box<dyn Error>> {
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
            primary_nodes.extend(p.vertices().iter().map(|&v| v.node_id));
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

    let used: HashSet<usize> = linear_blocks.iter().flat_map(|b| b.all_nodes()).collect();

    if let Some(output) = paths_fn {
        info!("Outputting paths in {}", output);
        let mut output = File::create(output)?;

        writeln!(output, "name\tlen\tpath\tassignment")?;

        for (block_id, block) in linear_blocks.into_iter().enumerate() {
            writeln!(
                output,
                "primary_{}\t{}\t{}\tPRIMARY",
                block_id,
                block.instance_path().total_length(&g),
                block.instance_path().print_format(&g, gaf_paths)
            )?;
            for (alt_id, &known_alt) in block.known_alt_nodes().iter().enumerate() {
                writeln!(
                    output,
                    "alt_{}_{}\t{}\t{}\tALT",
                    block_id,
                    alt_id,
                    g.node(known_alt).length,
                    Path::new(Vertex::forward(known_alt)).print_format(&g, gaf_paths)
                )?;
            }
        }

        for (node_id, n) in g.all_nodes().enumerate() {
            if !used.contains(&node_id) {
                writeln!(
                    output,
                    "unused_{}\t{}\t{}\tNA",
                    n.name,
                    n.length,
                    Path::new(Vertex::forward(node_id)).print_format(&g, gaf_paths)
                )?;
            }
        }
    }

    info!("All done");
    Ok(())
}

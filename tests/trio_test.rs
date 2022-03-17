extern crate log;

use graph_analysis::*;
use graph_analysis::trio::*;
use std::fs;

fn init() {
    let _ = env_logger::builder().is_test(true).try_init();
}

#[test]
fn homozygous_assignment() {
    init();

    let graph_fn = "tests/test_graphs/test1.gfa";
    let assignments_fn = "tests/test_graphs/test1.no_homozygous.csv";
    let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
    let assignments = trio::parse_read_assignments(&g, assignments_fn).unwrap();
    let assignments = trio::assign_homozygous(&g, assignments, 200_000,
        -1., usize::MAX);
    let mut homozygous_names : Vec<&str> = (0..g.node_cnt())
                                .filter(|&node_id| assignments.group(node_id) == Some(TrioGroup::HOMOZYGOUS))
                                .map(|node_id| g.name(node_id)).collect();
    homozygous_names.sort();
    assert_eq!(&homozygous_names, &["utig4-1237", "utig4-1552", "utig4-1826", "utig4-2589"]);
}

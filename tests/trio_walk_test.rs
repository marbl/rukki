#[macro_use]
extern crate log;

use graph_analysis::*;
use graph_analysis::graph::*;
use graph_analysis::trio::*;
use std::fs;

//#[test]
//fn test_folder() {
//    let path = env::current_dir().unwrap();
//    println!("The current directory is {}", path.display());
//    assert!(false);
//}

fn from_assignment_iterator<'a>(g: &'a Graph, node_assign_it: impl Iterator<Item=(usize, TrioGroup)>)
-> AssignmentStorage<'a> {
    let mut storage = AssignmentStorage::new(g);
    for (node_id, group) in node_assign_it {
        storage.update_group(node_id, group);
    }
    storage
}

fn from_parental_groups<'a>(g: &'a Graph, maternal: &[usize], paternal: &[usize])
-> AssignmentStorage<'a> {
    from_assignment_iterator(g, maternal.iter()
        .map(|n| (*n, TrioGroup::MATERNAL))
        .chain(paternal.iter()
            .map(|n| (*n, TrioGroup::PATERNAL))))
}

fn parse_group(group_str: &str) -> TrioGroup {
    match group_str {
        "MATERNAL" => TrioGroup::MATERNAL,
        "PATERNAL" => TrioGroup::PATERNAL,
        "HOMOZYGOUS" => TrioGroup::HOMOZYGOUS,
        _ => panic!("Invalid group string {group_str}"),
    }
}

fn parse_read_assignments<'a>(g: &'a Graph, assignments_fn: &str, load_homozygous: bool)
-> std::io::Result<trio::AssignmentStorage<'a>> {
    let mut assignments = trio::AssignmentStorage::new(g);
    for line in std::fs::read_to_string(assignments_fn)?.lines() {
        let split: Vec<&str> = line.trim().split('\t').collect();
        if &split[0].to_lowercase() != "node" && &split[0].to_lowercase() != "contig" {
            let node_name = split[0];
            let group = parse_group(split[1]);
            if load_homozygous || group != TrioGroup::HOMOZYGOUS {
                assignments.update_group(g.name2id(node_name), group);
            }
        }
    }
    Ok(assignments)
}

fn init() {
    let _ = env_logger::builder().is_test(true).try_init();
}

#[test]
fn assign_homozygous() {
    init();

    let graph_fn = "tests/test_graphs/test1.gfa";
    let assignments_fn = "tests/test_graphs/test1.ann";
    let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
    let assignments = parse_read_assignments(&g, assignments_fn, false).unwrap();
    let assignments = trio_walk::assign_homozygous(&g, assignments, 100_000);
    let ref_assignments = parse_read_assignments(&g, assignments_fn, true).unwrap();
    for (node_id, _) in g.all_nodes().enumerate() {
        assert!(assignments.group(node_id) == ref_assignments.group(node_id));
    }
}

#[test]
fn assign_homozygous_2() {
    init();

    let graph_fn = "tests/test_graphs/test1.gfa";
    let assignments_fn = "tests/test_graphs/test1.ann";
    let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
    let assignments = parse_read_assignments(&g, assignments_fn, false).unwrap();
    let assignments = trio_walk::assign_homozygous(&g, assignments, 50_000);
    let mut homozygous_names : Vec<&str> = (0..g.node_cnt())
                                .filter(|&node_id| assignments.group(node_id) == Some(TrioGroup::HOMOZYGOUS))
                                .map(|node_id| g.name(node_id)).collect();
    homozygous_names.sort();
    assert_eq!(&homozygous_names, &["utig4-1237", "utig4-1552", "utig4-1826", "utig4-2589"]);
}

#[test]
fn haplo_paths() {
    init();

    let graph_fn = "tests/test_graphs/test1.gfa";
    let assignments_fn = "tests/test_graphs/test1.ann";
    let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
    let assignments = parse_read_assignments(&g, assignments_fn, true).unwrap();

    let mut haplo_searcher = trio_walk::HaploSearcher::new(&g, &assignments, 500_000);
    let mut answer: Vec<(TrioGroup, String)> = haplo_searcher.find_all().into_iter()
                                           .map(|(p, _, group)| (group, p.print(&g))).collect();
    answer.sort();
    assert_eq!(&answer, &[
        (TrioGroup::MATERNAL,
            String::from("utig4-1829-,utig4-1826-,utig4-1828+,utig4-1832+,utig4-1245-,utig4-1240-,utig4-1237-,utig4-1239+,utig4-1552+,utig4-1554+,utig4-4105+,utig4-2593-,utig4-2589-,utig4-2590+")),
        (TrioGroup::PATERNAL,
            String::from("utig4-1830-,utig4-1826-,utig4-1827+,utig4-1831+,utig4-1243-,utig4-1241-,utig4-1237-,utig4-1238+,utig4-1552+,utig4-1553+,utig4-4096-,utig4-4097+,utig4-2592-,utig4-2589-,utig4-2591+"))])
}
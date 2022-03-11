extern crate log;

use graph_analysis::*;
use graph_analysis::trio::*;
use std::fs;

//fn from_assignment_iterator<'a>(g: &'a Graph, node_assign_it: impl Iterator<Item=(usize, TrioGroup)>)
//-> AssignmentStorage<'a> {
//    let mut storage = AssignmentStorage::new(g);
//    for (node_id, group) in node_assign_it {
//        storage.update_group(node_id, group);
//    }
//    storage
//}

//fn from_parental_groups<'a>(g: &'a Graph, maternal: &[usize], paternal: &[usize])
//-> AssignmentStorage<'a> {
//    from_assignment_iterator(g, maternal.iter()
//        .map(|n| (*n, TrioGroup::MATERNAL))
//        .chain(paternal.iter()
//            .map(|n| (*n, TrioGroup::PATERNAL))))
//}

fn init() {
    let _ = env_logger::builder().is_test(true).try_init();
}

#[test]
fn haplo_paths() {
    init();

    let graph_fn = "tests/test_graphs/test1.gfa";
    let assignments_fn = "tests/test_graphs/test1.ann.csv";
    let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
    let assignments = trio::parse_read_assignments(&g, assignments_fn).unwrap();

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


#[test]
fn augment_by_search() {
    init();

    let graph_fn = "tests/test_graphs/sparse_markers.gfa";
    let assignments_fn = "tests/test_graphs/sparse_markers.ann.csv";
    let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
    let assignments = trio::parse_read_assignments(&g, assignments_fn).unwrap();

    let init_node_len_thr = 500_000;
    assert_eq!(assignments.assigned().count(), 14);

    let augment_assign = augment_by_path_search(&g, assignments, init_node_len_thr);

    assert_eq!(augment_assign.assigned().count(), 17);
    assert_eq!(augment_assign.group(g.name2id("utig4-1421")), Some(TrioGroup::PATERNAL));
    assert_eq!(augment_assign.group(g.name2id("utig4-793")), Some(TrioGroup::PATERNAL));
    assert_eq!(augment_assign.group(g.name2id("utig4-1436")), Some(TrioGroup::PATERNAL));

    let mut haplo_searcher = trio_walk::HaploSearcher::new(&g,
        &augment_assign, init_node_len_thr);

    let mut answer: Vec<(TrioGroup, String)> = haplo_searcher.find_all().into_iter()
                                           .map(|(p, _, group)| (group, p.print(&g))).collect();
    answer.sort();
    assert_eq!(&answer, &[
        (TrioGroup::MATERNAL,
            String::from("utig4-1424-,utig4-1422-,utig4-1420-,utig4-1418+,utig4-792-,utig4-791+,utig4-795+,utig4-1435+,utig4-1437+,utig4-1439+")),
        (TrioGroup::PATERNAL,
            String::from("utig4-1423-,utig4-1421-,utig4-1419-,utig4-1418+,utig4-793-,utig4-791+,utig4-794+,utig4-1435+,utig4-1436+,utig4-1438+"))])
}

#[test]
fn bubble_filling() {
    init();

    let graph_fn = "tests/test_graphs/path_closing.gfa";
    let assignments_fn = "tests/test_graphs/path_closing.ann.csv";
    let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
    let assignments = trio::parse_read_assignments(&g, assignments_fn).unwrap();

    let init_node_len_thr = 500_000;
    assert_eq!(assignments.assigned().count(), 26);

    let augment_assign = augment_by_path_search(&g, assignments, init_node_len_thr);

    assert_eq!(augment_assign.assigned().count(), 28);
    assert_eq!(augment_assign.group(g.name2id("utig4-1397")), Some(TrioGroup::MATERNAL));
    assert_eq!(augment_assign.group(g.name2id("utig4-1347")), Some(TrioGroup::MATERNAL));

    let mut haplo_searcher = trio_walk::HaploSearcher::new(&g,
        &augment_assign, init_node_len_thr);

    haplo_searcher.try_fill_bubbles();

    let mut answer: Vec<(TrioGroup, String)> = haplo_searcher.find_all().into_iter()
                                           .map(|(p, _, group)| (group, p.print(&g))).collect();
    answer.sort();
    assert_eq!(&answer, &[
        (TrioGroup::MATERNAL,
            String::from("utig4-1575-,utig4-1574+,utig4-1397-,utig4-1395-,utig4-1347-,utig4-1343-,utig4-1344+,utig4-1568-,utig4-815-,utig4-814+,utig4-819+,utig4-1799-,utig4-1796-,utig4-1798+")),
        (TrioGroup::MATERNAL,
            String::from("utig4-3444+,utig4-4080-,utig4-771-,utig4-768-,utig4-770+")),
        (TrioGroup::PATERNAL,
            String::from("utig4-1576-,utig4-1574+,utig4-1396-,utig4-1395-,utig4-1346-,utig4-1343-,utig4-1344+,utig4-1568-,utig4-815-,utig4-814+,utig4-818+,utig4-1796-,utig4-1797+")),
        (TrioGroup::PATERNAL,
            String::from("utig4-3412+,utig4-774-,utig4-772-,utig4-768-,utig4-769+"))]);
}

#[test]
fn haplo_paths_2() {
    init();

    let graph_fn = "tests/test_graphs/test2.gfa";
    let assignments_fn = "tests/test_graphs/test2.ann.csv";
    let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
    let assignments = trio::parse_read_assignments(&g, assignments_fn).unwrap();

    let init_node_len_thr = 500_000;
    assert_eq!(assignments.assigned().count(), 42);

    let augment_assign = augment_by_path_search(&g, assignments, init_node_len_thr);

    assert_eq!(augment_assign.group(g.name2id("utig4-414")), Some(TrioGroup::MATERNAL));
    assert_eq!(augment_assign.group(g.name2id("utig4-308")), Some(TrioGroup::MATERNAL));
    assert_eq!(augment_assign.group(g.name2id("utig4-415")), Some(TrioGroup::PATERNAL));
    assert_eq!(augment_assign.assigned().count(), 45);

    let mut haplo_searcher = trio_walk::HaploSearcher::new(&g,
        &augment_assign, init_node_len_thr);

    let mut answer: Vec<(TrioGroup, String)> = haplo_searcher.find_all().into_iter()
                                           .map(|(p, _, group)| (group, p.print(&g))).collect();
    answer.sort();
    assert_eq!(&answer, &[
        (TrioGroup::MATERNAL,
            String::from("utig4-3444+,utig4-4080-,utig4-771-,utig4-768-,utig4-770+,utig4-1384-,utig4-1385+,utig4-1898-,utig4-1897-,utig4-414-,utig4-412+,utig4-416+,utig4-419+,utig4-4073-,utig4-1460-,utig4-1459+,utig4-1463+,utig4-4227+,utig4-311-,utig4-307-,utig4-308+,utig4-3431-,utig4-3430+")),
        (TrioGroup::PATERNAL,
            String::from("utig4-3412+,utig4-774-,utig4-772-,utig4-768-,utig4-769+,utig4-1384-,utig4-1386+,utig4-1899-,utig4-1897-,utig4-413-,utig4-412+,utig4-415+,utig4-422-,utig4-418+,utig4-421+,utig4-1461-,utig4-1459+,utig4-1462+,utig4-4227+,utig4-312-,utig4-307-,utig4-309+,utig4-3429+"))]);
}

#[test]
fn haplo_paths_3() {
    init();

    let graph_fn = "tests/test_graphs/test3.gfa";
    let assignments_fn = "tests/test_graphs/test3.ann.csv";
    let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
    let assignments = trio::parse_read_assignments(&g, assignments_fn).unwrap();

    let init_node_len_thr = 500_000;
    assert_eq!(assignments.assigned().count(), 76);

    let augment_assign = augment_by_path_search(&g, assignments, init_node_len_thr);

    assert_eq!(augment_assign.group(g.name2id("utig4-1404")), Some(TrioGroup::PATERNAL));
    assert_eq!(augment_assign.group(g.name2id("utig4-1403")), Some(TrioGroup::MATERNAL));

    assert_eq!(augment_assign.assigned().count(), 82);

    let mut haplo_searcher = trio_walk::HaploSearcher::new(&g,
        &augment_assign, init_node_len_thr);

    let mut answer: Vec<(TrioGroup, String)> = haplo_searcher.find_all().into_iter()
                                           .map(|(p, _, group)| (group, p.print(&g))).collect();
    answer.sort();
    assert_eq!(&answer, &[
        (TrioGroup::MATERNAL,
            String::from("utig4-4093-,utig4-3587-,utig4-3588+,utig4-4041-,utig4-3592+,utig4-1535-,utig4-1533-,utig4-1529-,utig4-1531+,utig4-1892-,utig4-925-,utig4-923+,utig4-926+,utig4-1595+,utig4-1597+,utig4-1896+,utig4-1619-,utig4-1617+,utig4-65-,utig4-64+,utig4-67+,AMBIG,utig4-1477-,utig4-1476+,AMBIG,utig4-1249+,utig4-1252+,utig4-1254+,utig4-3626+,utig4-3631+,utig4-1027-,utig4-1025-,utig4-1022-,utig4-1019-,utig4-1020+,utig4-1387+,AMBIG,utig4-1392+,utig4-1393+,utig4-1450+,utig4-1451+,utig4-1795+,utig4-1406-,utig4-1402-,utig4-1403+,utig4-3448-,utig4-1409+,utig4-3446-,GAP,utig4-3456-")),
        (TrioGroup::PATERNAL,
            String::from("utig4-3455-,utig4-3445-,utig4-3447+,utig4-1410-,utig4-1408-,utig4-1404-,utig4-1402+,utig4-1405+,utig4-1795-,utig4-1452-,utig4-1450-,utig4-1394-,utig4-1392-,AMBIG,utig4-1387-,utig4-1021-,utig4-1019+,utig4-1023+,utig4-1024+,utig4-1026+,utig4-3630-,utig4-3626-,utig4-3627+,utig4-1257-,utig4-1253-,utig4-1249-,AMBIG,utig4-1476-,utig4-1478+,utig4-3650-,utig4-68-,utig4-64-,utig4-66+,utig4-1617-,utig4-1618+,utig4-1896-,utig4-1596-,utig4-1595-,utig4-927-,utig4-923-,utig4-924+,utig4-1892+,utig4-1530-,utig4-1529+,utig4-1532+,utig4-1534+,utig4-3593-,utig4-3591-,utig4-3589-,GAP,utig4-3384+"))]);
}
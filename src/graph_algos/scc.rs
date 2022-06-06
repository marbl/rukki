use crate::graph::*;
use super::dfs;
use super::only_or_none;
use std::collections::{HashSet,HashMap};
use log::debug;

//Implementing Kosaraju-Sharir algorithm
//'trivial' SCCs of individual vertices are not reported
//NB. Loop of single vertex is considered 'NON-trivial'
pub fn strongly_connected(graph: &Graph) -> Vec<Vec<Vertex>> {
    let mut non_trivial_sccs: Vec<Vec<Vertex>> = Vec::new();
    let is_loop = |v: Vertex| {
        graph.outgoing_edges(v).iter().any(|l| l.end == v)
    };

    // run DFS on direct edges
    let mut dfs = dfs::DFS::new_forward(graph);
    dfs.run();
    let mut used: HashSet<Vertex> = HashSet::new();
    // consider vertices in decreasing order of exit times (latest exit times first)
    for &v in dfs.exit_order().iter().rev() {
        if !used.contains(&v) {
            // run DFS on reverse edges
            let mut reverse_dfs = dfs::DFS::new_reverse(graph);
            reverse_dfs.set_blocked(used);
            reverse_dfs.run_from(v);
            let visited = reverse_dfs.exit_order();
            assert!(!visited.is_empty());
            if visited.len() > 1 || is_loop(visited[0]) {
                debug!("Identified non-trivial component of size {}: {}",
                        visited.len(),
                        visited.iter().map(|&v| graph.v_str(v)).collect::<Vec<String>>().join(","));

                non_trivial_sccs.push(visited.clone());
            }
            used = reverse_dfs.take_blocked();
        }
    }
    assert!(check_consistency(graph, &non_trivial_sccs));
    non_trivial_sccs
}

pub fn nodes_in_sccs(_g: &Graph, sccs: &[Vec<Vertex>]) -> HashSet<usize> {
    HashSet::from_iter(sccs.iter()
        .flat_map(|comp| comp.iter().map(|v| v.node_id)))
}

fn check_consistency(graph: &Graph, non_trivial_sccs: &[Vec<Vertex>]) -> bool {
    let mut vertices_to_scc = HashMap::new();
    for (scc_id, vertices) in non_trivial_sccs.iter().enumerate() {
        for v in vertices {
            vertices_to_scc.insert(v, scc_id);
        }
    }

    let mut considered_node_ids : HashSet<usize> = HashSet::new();
    for v in graph.all_vertices() {
        if considered_node_ids.contains(&v.node_id) {
            continue;
        }
        let sorted = |mut vertices: Vec<Vertex>| {
            vertices.sort();
            vertices
        };
        if let Some(&scc_id) = vertices_to_scc.get(&v) {
            for scc_v in &non_trivial_sccs[scc_id] {
                considered_node_ids.insert(scc_v.node_id);
            }
            match vertices_to_scc.get(&v.rc()) {
                None => return false,
                Some(&rc_scc_id) => {
                    assert_eq!(sorted(non_trivial_sccs[scc_id].clone()),
                        sorted(non_trivial_sccs[rc_scc_id].iter().map(|w| w.rc()).collect()));
                },
            }
        }
    }
    true
}

//Building condensation Graph
pub fn condensation(graph: &Graph, non_trivial_sccs: &[Vec<Vertex>], ignore_loops: bool) -> (Graph, HashMap<Vertex, Vertex>) {
    assert!(check_consistency(graph, non_trivial_sccs));
    let mut condensation = Graph::new();
    let mut vertices_to_scc = HashMap::new();
    for (scc_id, vertices) in non_trivial_sccs.iter().enumerate() {
        //filtering 'trivial' loops
        if vertices.len() == 1 {
            continue;
        }
        for v in vertices {
            vertices_to_scc.insert(v, scc_id);
        }
    }

    let mut old_2_new : HashMap<Vertex, Vertex> = HashMap::new();

    let mut update_old_2_new = |old_vertices: &[Vertex], new_node_id: usize| {
        //two passes for more consistent processing of self-conjugate scc
        for v in old_vertices {
            old_2_new.insert(v.rc(), Vertex::reverse(new_node_id));
        }
        for v in old_vertices {
            old_2_new.insert(*v, Vertex::forward(new_node_id));
        }
    };

    let mut considered_node_ids : HashSet<usize> = HashSet::new();
    for (node_id, node) in graph.node_iter().enumerate() {
        let v = Vertex::forward(node_id);
        if considered_node_ids.contains(&node_id) {
            continue;
        }
        if let Some(&scc_id) = vertices_to_scc.get(&v) {
            let scc_vertices = &non_trivial_sccs[scc_id];
            for scc_v in scc_vertices {
                considered_node_ids.insert(scc_v.node_id);
            }
            if scc_vertices.contains(&v.rc()) {
                debug!("Dealing with self-conjugate SCC {}: {}", scc_id,
                    scc_vertices.iter().map(|&w| graph.v_str(w)).collect::<Vec<String>>().join(""))
            }
            let length = scc_vertices.iter().map(|w| graph.node(w.node_id).length).max().unwrap();
            let name = format!("scc_{}_vcnt_{}_init_{}", scc_id, scc_vertices.len(), node.name);
            //let cnd_node;
            let cnd_id = condensation.add_node(Node{name, length, coverage: 0.,});
            update_old_2_new(scc_vertices, cnd_id);
        } else {
            considered_node_ids.insert(v.node_id);
            let cnd_id = condensation.add_node(node.clone());
            update_old_2_new(std::slice::from_ref(&v), cnd_id);
        }
    }

    for l in graph.all_links() {
        let &v = old_2_new.get(&l.start).unwrap();
        let &w = old_2_new.get(&l.end).unwrap();
        //checking that no link between nodes exists
        if ignore_loops && v == w {
            debug!("Loop ignored for vertex {}", condensation.v_str(v));
            continue
        }
        if !condensation.outgoing_edges(v).iter().any(|l| l.end == w) {
            condensation.add_link(Link {
                start: v,
                end: w,
                overlap: l.overlap,
            });
        }
    }

    (condensation, old_2_new)
}

pub struct LocalizedTangle {
    pub entrance: Link,
    pub exit: Link,
    pub vertices: Vec<Vertex>,
}

//very crude (under-)estimate without multiplicity guessing!
//subtracts minimal incoming overlap from every vertex and takes sum
pub fn estimate_size_no_mult(tangle: &LocalizedTangle, g: &Graph) -> usize {
    let shortest_incoming_overlap = |v: Vertex| {
        return g.incoming_edges(v).iter().map(|l| l.overlap).min().unwrap_or(0);
    };

    tangle.vertices.iter()
        .map(|&v| g.vertex_length(v) - shortest_incoming_overlap(v))
        .sum()
}

fn find_localized(g: &Graph, non_trivial_scc: &[Vertex]) -> Option<LocalizedTangle> {
    let component_vertices: HashSet<Vertex>
        = HashSet::from_iter(non_trivial_scc.iter().copied());

    //TODO learn if there is a better way or implement my own helper function
    let entrance = only_or_none(component_vertices.iter()
                        .flat_map(|&v| g.incoming_edges(v))
                        .filter(|l| !component_vertices.contains(&l.start)))?;

    let exit = only_or_none(component_vertices.iter()
                        .flat_map(|&v| g.outgoing_edges(v))
                        .filter(|l| !component_vertices.contains(&l.end)))?;

    //TODO think where this check should be performed
    //Also checking that entrance and exit
    //are the only ways to go from corresponding vertices
    let entrance = only_or_none(g.outgoing_edges(entrance.start).into_iter())?;
    let exit = only_or_none(g.incoming_edges(exit.end).into_iter())?;

    //guard against potential tricky strand-switching case
    if entrance.start.node_id == exit.end.node_id {
        None
    } else {
        Some(LocalizedTangle {
            entrance,
            exit,
            vertices : non_trivial_scc.to_owned(),
        })
    }
}

pub fn find_small_localized(g: &Graph,
    non_trivial_sccs: &[Vec<Vertex>],
    size_limit: usize)
-> Vec<LocalizedTangle> {
    non_trivial_sccs.iter()
        .filter_map(|vs| find_localized(g, vs))
        .filter(|t| estimate_size_no_mult(t, g) <= size_limit)
        .collect()
}
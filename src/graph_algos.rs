pub mod scc;
pub mod superbubble;

pub mod graph_algos {
use crate::graph::*;
use std::collections::HashSet;

fn inner_dfs(g: &Graph, v: Vertex, node_len_thr: usize, visited: &mut HashSet<Vertex>, long_ext: &mut Vec<Vertex>) {
    visited.insert(v);
    //if only one vertex is visited then it means we just started
    if visited.len() > 1 && g.node(v.node_id).length >= node_len_thr {
        long_ext.push(v);
    } else {
        for l in g.outgoing_edges(v) {
            let w = l.end;
            if !visited.contains(&w) {
                inner_dfs(g, w, node_len_thr, visited, long_ext);
            }
        }
    }
}

pub fn bounded_dfs(g: &Graph, v: Vertex, node_len_thr: usize)
    -> (HashSet<Vertex>, Vec<Vertex>) {
    //TODO change for integer vectors
    let mut visited = HashSet::new();
    let mut long_ext = Vec::new();
    inner_dfs(g, v, node_len_thr, &mut visited, &mut long_ext);
    (visited, long_ext)
}
}

use crate::graph::*;
use std::collections::HashSet;
use log::debug;

#[derive(Copy, Clone, Debug)]
pub enum TraversalDirection {
    FORWARD,
    REVERSE,
    //TODO add
    //UNDIRECTED,
}

//TODO replace sets with arrays

//FIXME pass functions giving neighbour iterators to allow more flexibility (directions, subgraphs, length boundaries, etc)
//TODO use within trio_walk
struct DFS<'a> {
    g: &'a Graph,
    direction: TraversalDirection,
    blocked: HashSet<Vertex>,
    tout: Vec<Vertex>,
}

impl<'a> DFS<'a> {

    pub fn new(g: &'a Graph, direction: TraversalDirection, blocked: HashSet<Vertex>) -> DFS<'a> {
        DFS {
            g,
            direction,
            blocked,
            tout: Vec::new(),
        }
    }

    pub fn new_forward(g: &'a Graph) -> DFS<'a> {
        Self::new(g, TraversalDirection::FORWARD, HashSet::new())
    }

    pub fn new_reverse(g: &'a Graph) -> DFS<'a> {
        Self::new(g, TraversalDirection::REVERSE, HashSet::new())
    }

    //TODO use iterators
    fn neighbors(&self, v: Vertex) -> Vec<Vertex> {
        match self.direction {
            TraversalDirection::FORWARD => self.g.outgoing_edges(v).iter()
                                                 .map(|l| l.end).collect(),
            TraversalDirection::REVERSE => self.g.incoming_edges(v).iter()
                                                 .map(|l| l.start).collect(),
        }
    }

    pub fn run_from(&mut self, v: Vertex) {
        self.blocked.insert(v);

        for w in self.neighbors(v) {
            if !self.blocked.contains(&w) {
                self.run_from(w);
            }
        }

        self.tout.push(v);
    }

    pub fn run(&mut self) {
        for v in self.g.all_vertices() {
            if !self.blocked.contains(&v) {
                self.run_from(v);
            }
        }
    }

    pub fn blocked(self) -> HashSet<Vertex> {
        self.blocked
    }

    pub fn exit_order(&self) -> &Vec<Vertex> {
        &self.tout
    }
}

//Implementing Kosaraju-Sharir algorithm
//'trivial' SCCs are not reported (loop of single vertex is not 'trivial')
pub fn strongly_connected(graph: &Graph) -> Vec<Vec<Vertex>> {
    let mut non_trivial_sccs: Vec<Vec<Vertex>> = Vec::new();
    let is_loop = |v: Vertex| {
        graph.outgoing_edges(v).iter().any(|l| l.end == v)
    };

    // run DFS on direct edges
    let mut dfs = DFS::new_forward(graph);
    dfs.run();
    let mut used: HashSet<Vertex> = HashSet::new();
    // consider vertices in decreasing order of exit times (latest exit times first)
    for &v in dfs.exit_order().iter().rev() {
        if !used.contains(&v) {
            // run DFS on reverse edges
            let mut reverse_dfs = DFS::new(graph, TraversalDirection::REVERSE, used);
            reverse_dfs.run_from(v);
            let reached = reverse_dfs.exit_order();
            assert!(reached.len() > 0);
            if reached.len() > 1 || is_loop(reached[0]) {
                debug!("Identified non-trivial component of size {}: {}",
                        reached.len(),
                        reached.iter().map(|&v| graph.v_str(v)).collect::<Vec<String>>().join(","));

                non_trivial_sccs.push(reached.clone());
            }
            used = reverse_dfs.blocked();
        }
    }
    non_trivial_sccs
}
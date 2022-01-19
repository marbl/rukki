use crate::graph::*;
use std::collections::{HashSet,HashMap};
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
    assert!(check_consistency(graph, &non_trivial_sccs));
    non_trivial_sccs
}

fn check_consistency(graph: &Graph, non_trivial_sccs: &Vec<Vec<Vertex>>) -> bool {
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
pub fn condensation(graph: &Graph, non_trivial_sccs: &Vec<Vec<Vertex>>, ignore_loops: bool) -> (Graph, HashMap<Vertex, Vertex>) {
    assert!(check_consistency(graph, &non_trivial_sccs));
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
            update_old_2_new(&scc_vertices, cnd_id);
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
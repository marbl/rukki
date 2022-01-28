use crate::graph::*;
use std::collections::HashSet;

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
pub struct DFS<'a> {
    g: &'a Graph,
    direction: TraversalDirection,
    blocked: HashSet<Vertex>,
    tout: Vec<Vertex>,
    node_len_thr: usize,
}

impl<'a> DFS<'a> {

    pub fn new(g: &'a Graph, direction: TraversalDirection) -> DFS<'a> {
        DFS {
            g,
            direction,
            blocked: HashSet::new(),
            tout: Vec::new(),
            node_len_thr: usize::MAX,
        }
    }

    pub fn new_forward(g: &'a Graph) -> DFS<'a> {
        Self::new(g, TraversalDirection::FORWARD)
    }

    pub fn new_reverse(g: &'a Graph) -> DFS<'a> {
        Self::new(g, TraversalDirection::REVERSE)
    }

    //FIXME make consume self and return new DFS
    pub fn set_blocked(&mut self, blocked: HashSet<Vertex>) {
        self.blocked = blocked;
    }

    //FIXME make consume self and return new DFS
    pub fn set_max_node_len(&mut self, max_node_len: usize) {
        self.node_len_thr = max_node_len;
    }

    //FIXME make consume self and return new DFS
    pub fn extend_blocked(&mut self, iter: impl IntoIterator<Item=Vertex>) {
        self.blocked.extend(iter);
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
        assert!(!self.blocked.contains(&v));
        self.blocked.insert(v);

        for w in self.neighbors(v) {
            if !self.blocked.contains(&w)
                && self.g.vertex_length(w) <= self.node_len_thr {
                self.run_from(w);
            }
        }

        self.tout.push(v);
    }

    //todo maybe rename into topsort?
    //will run from long nodes, but not from blocked
    pub fn run(&mut self) {
        assert!(self.node_len_thr == usize::MAX);
        for v in self.g.all_vertices() {
            if !self.blocked.contains(&v) {
                self.run_from(v);
            }
        }
    }

    //includes visited and initially blocked
    pub fn take_blocked(self) -> HashSet<Vertex> {
        self.blocked
    }

    pub fn blocked(&self) -> &HashSet<Vertex> {
        &self.blocked
    }

    pub fn exit_order(&self) -> &Vec<Vertex> {
        &self.tout
    }

    //todo return iterator
    //nodes that were reached, but not visited (initally blocked or too long)
    pub fn boundary(&self) -> Vec<Vertex> {
        let mut boundary = Vec::new();
        let visited : HashSet<Vertex> = self.tout.iter().copied().collect();

        for &v in &visited {
            for w in self.neighbors(v) {
                if !visited.contains(&w) {
                    boundary.push(w);
                }
            }
        }
        boundary
    }

}

//includes boundary vertices (longer than threshold) and visited dead-ends
//currently will include v if it's a dead-end (but not if it's longer than threshold)
pub fn sinks_ahead(g: &Graph, v: Vertex, node_len_thr: usize) -> Vec<Vertex> {
    let mut dfs = DFS::new_forward(g);
    dfs.set_max_node_len(node_len_thr);
    //inner_dfs(g, v, node_len_thr, &mut visited, &mut border);
    dfs.run_from(v);
    let mut sinks = dfs.boundary();
    assert!(sinks.iter().all(|&x| g.vertex_length(x) >= node_len_thr));
    sinks.extend(dfs.exit_order().iter().filter(|&x| g.outgoing_edge_cnt(*x) == 0).copied());
    sinks
}
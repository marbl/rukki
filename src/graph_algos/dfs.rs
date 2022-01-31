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

struct ShortNodeComponent {
    sources: HashSet<Vertex>,
    sinks: HashSet<Vertex>,
    has_deadends: bool,
    reached: HashSet<Vertex>,
}

impl ShortNodeComponent {

    fn consider(&mut self, g: &Graph, v: Vertex, l: Link, length_threshold: usize) {
        if !self.reached.contains(&v) {
            self.reached.insert(v);

            if g.vertex_length(v) >= length_threshold
                && v == l.start {
                //if v is long and we came into its end then it's a source
                self.sources.insert(v);
            } else {
                //otherwise consider it's incoming edges
                if g.incoming_edge_cnt(v) == 0 {
                    self.has_deadends = true;
                }
                for i_l in g.incoming_edges(v) {
                    if i_l != l {
                        self.consider(g, i_l.start, i_l, length_threshold);
                    }
                }
            }

            if g.vertex_length(v) >= length_threshold
                && v == l.end {
                //if v is long and we came into its start then it's a sink
                self.sinks.insert(v);
            } else {
                //otherwise consider it's outgoing edges
                if g.outgoing_edge_cnt(v) == 0 {
                    self.has_deadends = true;
                }
                for o_l in g.outgoing_edges(v) {
                    if o_l != l {
                        self.consider(g, o_l.end, o_l, length_threshold);
                    }
                }
            }
        }
    }

    //returns true if all nodes are distinct within sources/sinks union
    fn simple_boundary(&self) -> bool {
        let mut used = HashSet::new();
        for v in self.sinks.iter().chain(self.sources.iter()) {
            if used.contains(&v.node_id) {
                return false;
            }
            used.insert(v.node_id);
        }
        return true;
    }

    pub fn ahead_from_long(g: &Graph, v: Vertex, length_threshold: usize) -> ShortNodeComponent {
        assert!(g.vertex_length(v) >= length_threshold);
        let mut component = ShortNodeComponent {
            sources: std::iter::once(v).collect(),
            sinks: HashSet::new(),
            has_deadends: (g.outgoing_edge_cnt(v) == 0),
            reached: std::iter::once(v).collect(),
        };

        for o_l in g.outgoing_edges(v) {
            component.consider(g, o_l.end, o_l, length_threshold);
        }
        component
    }
}
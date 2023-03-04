use crate::graph::*;
use itertools::Itertools;
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
    consider_f: Option<&'a dyn Fn(Vertex) -> bool>,
    blocked: HashSet<Vertex>,
    tout: Vec<Vertex>,
    node_len_thr: usize,
}

impl<'a> DFS<'a> {
    pub fn new(
        g: &'a Graph,
        direction: TraversalDirection,
        consider_f: Option<&'a dyn Fn(Vertex) -> bool>,
    ) -> DFS<'a> {
        DFS {
            g,
            direction,
            consider_f,
            blocked: HashSet::new(),
            tout: Vec::new(),
            node_len_thr: usize::MAX,
        }
    }

    pub fn new_forward(g: &'a Graph) -> DFS<'a> {
        Self::new(g, TraversalDirection::FORWARD, None)
    }

    pub fn new_reverse(g: &'a Graph) -> DFS<'a> {
        Self::new(g, TraversalDirection::REVERSE, None)
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
    pub fn extend_blocked(&mut self, iter: impl IntoIterator<Item = Vertex>) {
        self.blocked.extend(iter);
    }

    //TODO use iterators
    fn neighbors(&self, v: Vertex) -> Vec<Vertex> {
        match self.direction {
            TraversalDirection::FORWARD => self
                .g
                .outgoing_edges(v)
                .iter()
                .map(|l| l.end)
                .filter(|&x| self.consider_f.is_none() || self.consider_f.unwrap()(x))
                .collect(),
            TraversalDirection::REVERSE => self
                .g
                .incoming_edges(v)
                .iter()
                .map(|l| l.start)
                .filter(|&x| self.consider_f.is_none() || self.consider_f.unwrap()(x))
                .collect(),
        }
    }

    pub fn run_from(&mut self, v: Vertex) {
        assert!(!self.blocked.contains(&v));
        self.blocked.insert(v);

        for w in self.neighbors(v) {
            if !self.blocked.contains(&w) && self.g.vertex_length(w) <= self.node_len_thr {
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
        let mut boundary = HashSet::new();
        let visited: HashSet<Vertex> = self.tout.iter().copied().collect();

        for &v in &visited {
            for w in self.neighbors(v) {
                if !visited.contains(&w) {
                    boundary.insert(w);
                }
            }
        }
        boundary.into_iter().collect()
    }

    //TODO return iterator?
    //return nodes that didn't have any neighbors
    pub fn dead_ends(&self) -> Vec<Vertex> {
        self.exit_order()
            .iter()
            .filter(|&v| self.neighbors(*v).is_empty())
            .copied()
            .collect()
    }
}

//includes boundary vertices (longer than threshold) and visited dead-ends
//currently will include v if it's a dead-end (but not if it's longer than threshold)
//returns pair of sinks and all ('inner') visited vertices
//visited vertices will overlap sinks by short dead-ends
pub fn sinks_ahead(
    g: &Graph,
    v: Vertex,
    node_len_thr: usize,
    subgraph_f: Option<&dyn Fn(Vertex) -> bool>,
) -> (Vec<Vertex>, HashSet<Vertex>) {
    let mut dfs = DFS::new(g, TraversalDirection::FORWARD, subgraph_f);
    dfs.set_max_node_len(node_len_thr);
    //inner_dfs(g, v, node_len_thr, &mut visited, &mut border);
    dfs.run_from(v);
    let mut sinks = dfs.boundary();
    assert!(sinks.iter().all(|&x| g.vertex_length(x) >= node_len_thr));
    //extend to dead-ends
    sinks.extend(dfs.dead_ends());
    (sinks, dfs.take_blocked())
}

//includes boundary vertices (longer than threshold) and visited dead-ends
//currently will include v if it's a dead-end (but not if it's longer than threshold)
//returns pair of sources and all ('inner') visited vertices
//visited vertices will overlap sources by short dead-ends
pub fn sources_behind(
    g: &Graph,
    v: Vertex,
    node_len_thr: usize,
    subgraph_f: Option<&dyn Fn(Vertex) -> bool>,
) -> (Vec<Vertex>, HashSet<Vertex>) {
    let mut dfs = DFS::new(g, TraversalDirection::REVERSE, subgraph_f);
    dfs.set_max_node_len(node_len_thr);
    //inner_dfs(g, v, node_len_thr, &mut visited, &mut border);
    dfs.run_from(v);
    let mut sources = dfs.boundary();
    assert!(sources.iter().all(|&x| g.vertex_length(x) >= node_len_thr));
    //extend to dead-ends
    sources.extend(dfs.dead_ends());
    (sources, dfs.take_blocked())
}

pub struct ShortNodeComponent {
    pub sources: HashSet<Vertex>,
    pub sinks: HashSet<Vertex>,
    pub has_deadends: bool,
    pub inner: HashSet<Vertex>,
}

impl ShortNodeComponent {
    fn consider(&mut self, g: &Graph, v: Vertex, l: Link, length_threshold: usize) {
        let mut is_source = false;
        let mut is_sink = false;

        if g.vertex_length(v) < length_threshold {
            if !self.inner.insert(v) {
                //inner already considered
                return;
            }
        } else {
            //v is long
            if v == l.start {
                //if v is long and we came from the 'right'
                is_source = true;
                if !self.sources.insert(v) {
                    //source already considered
                    return;
                }
            } else {
                assert!(v == l.end);
                //if v is long and we came from the 'left'
                is_sink = true;
                if !self.sinks.insert(v) {
                    //sink already considered
                    return;
                }
            }
        }

        //if not a source consider it's incoming edges
        if !is_source {
            if g.incoming_edge_cnt(v) == 0 {
                assert!(g.vertex_length(v) < length_threshold);
                self.has_deadends = true;
            }
            for i_l in g.incoming_edges(v) {
                if i_l != l {
                    self.consider(g, i_l.start, i_l, length_threshold);
                }
            }
        }

        //if not a sink consider outgoing edges
        if !is_sink {
            if g.outgoing_edge_cnt(v) == 0 {
                assert!(g.vertex_length(v) < length_threshold);
                self.has_deadends = true;
            }
            for o_l in g.outgoing_edges(v) {
                if o_l != l {
                    self.consider(g, o_l.end, o_l, length_threshold);
                }
            }
        }
    }

    //returns true if all nodes are distinct within sources/sinks union
    pub fn simple_boundary(&self) -> bool {
        let mut used = HashSet::new();
        for v in self.sinks.iter().chain(self.sources.iter()) {
            if used.contains(&v.node_id) {
                return false;
            }
            used.insert(v.node_id);
        }
        true
    }

    pub fn ahead_from_long(g: &Graph, v: Vertex, length_threshold: usize) -> ShortNodeComponent {
        assert!(g.vertex_length(v) >= length_threshold);
        let mut component = ShortNodeComponent {
            sources: std::iter::once(v).collect(),
            sinks: HashSet::new(),
            has_deadends: false,
            inner: HashSet::new(),
        };

        for o_l in g.outgoing_edges(v) {
            component.consider(g, o_l.end, o_l, length_threshold);
        }
        component
    }

    pub fn back_from_long(g: &Graph, v: Vertex, length_threshold: usize) -> ShortNodeComponent {
        assert!(g.vertex_length(v) >= length_threshold);
        let mut component = ShortNodeComponent {
            sources: HashSet::new(),
            sinks: std::iter::once(v).collect(),
            has_deadends: false,
            inner: HashSet::new(),
        };

        for i_l in g.incoming_edges(v) {
            component.consider(g, i_l.start, i_l, length_threshold);
        }
        component
    }

    //todo refactor and simplify logic!
    //if v is long searching ahead from it, otherwise search in both directions
    pub fn search_from(g: &Graph, v: Vertex, length_threshold: usize) -> ShortNodeComponent {
        if g.vertex_length(v) >= length_threshold {
            Self::ahead_from_long(g, v, length_threshold)
        } else {
            let mut component = ShortNodeComponent {
                sources: HashSet::new(),
                sinks: HashSet::new(),
                has_deadends: (g.outgoing_edge_cnt(v) == 0 || g.incoming_edge_cnt(v) == 0),
                inner: std::iter::once(v).collect(),
            };
            for i_l in g.incoming_edges(v) {
                component.consider(g, i_l.start, i_l, length_threshold);
            }
            for o_l in g.outgoing_edges(v) {
                component.consider(g, o_l.end, o_l, length_threshold);
            }
            component
        }
    }

    pub fn all_nodes(&self) -> impl Iterator<Item = &Vertex> {
        self.inner
            .iter()
            .chain(self.sources.iter())
            .chain(self.sinks.iter())
    }

    pub fn print(&self, g: &Graph) -> String {
        format!(
            "Sources: {}; sinks: {}",
            self.sources.iter().map(|&v| g.v_str(v)).join(", "),
            self.sinks.iter().map(|&v| g.v_str(v)).join(", ")
        )
    }
}

use crate::graph::*;
use std::cmp;
use std::collections::HashSet;
use std::collections::HashMap;
use log::debug;

type DistRange = (usize, usize);

fn shift_range((min, max): DistRange, s: usize) -> DistRange {
    (min + s, max + s)
}

fn merge_range((min1, max1): DistRange, (min2, max2): DistRange) -> DistRange {
    (cmp::min(min1, min2), cmp::max(max1, max2))
}

pub struct Superbubble {
    start_vertex: Vertex,
    end_vertex: Option<Vertex>,
    //vertex to path length range
    reached_vertices: HashMap<Vertex, DistRange>,
}

impl Superbubble {

    fn link_dist_range(&self, l: Link, g: &Graph) -> Option<DistRange> {
        let &r = self.reached_vertices.get(&l.start)?;
        let enode_len = g.vertex_length(l.end);
        assert!(enode_len >= l.overlap);
        Some(shift_range(r, enode_len - l.overlap))
    }

    pub fn longest_path(&self, g: &Graph) -> Path {
        let mut v = self.end_vertex.unwrap();
        let mut longest_dist = self.reached_vertices.get(&v).unwrap().1;
        let mut rc_p = Path::new(v.rc());
        'outer: while v != self.start_vertex {
            //let l = self.heaviest_backtrace.get(v).unwrap();
            for l in g.incoming_edges(v) {
                if let Some((_, l_d)) = self.link_dist_range(l, g) {
                    if l_d == longest_dist {
                        assert!(l.end == v);
                        rc_p.append(l.rc());
                        v = l.start;
                        longest_dist = self.reached_vertices.get(&l.start).unwrap().1;
                        continue 'outer;
                    }
                }
            }
            panic!("Couldn't recover bubble path");
        }
        rc_p.reverse_complement()
    }

    pub fn shortest_path(&self, g: &Graph) -> Path {
        let mut v = self.end_vertex.unwrap();
        let mut shortest_dist = self.reached_vertices.get(&v).unwrap().0;
        let mut rc_p = Path::new(v.rc());
        'outer: while v != self.start_vertex {
            //let l = self.heaviest_backtrace.get(v).unwrap();
            for l in g.incoming_edges(v) {
                if let Some((l_d, _)) = self.link_dist_range(l, g) {
                    if l_d == shortest_dist {
                        assert!(l.end == v);
                        rc_p.append(l.rc());
                        v = l.start;
                        shortest_dist = self.reached_vertices.get(&l.start).unwrap().0;
                        continue 'outer;
                    }
                }
            }
            panic!("Couldn't recover bubble path");
        }
        rc_p.reverse_complement()
    }

    pub fn vertices(&self) -> impl Iterator<Item=&Vertex> + '_ {
        self.reached_vertices.keys()
    }

    pub fn inner_vertices(&self) -> impl Iterator<Item=&Vertex> + '_ {
        self.reached_vertices.keys().filter(|&v| *v != self.start_vertex() && *v != self.end_vertex())
    }

    pub fn start_vertex(&self) -> Vertex {
        self.start_vertex
    }

    pub fn end_vertex(&self) -> Vertex {
        self.end_vertex.unwrap()
    }

    pub fn length_range(&self, g: &Graph) -> (usize, usize) {
        let r = *self.reached_vertices.get(&self.end_vertex()).unwrap();
        //currently start vertex and end vertex can't be the same
        assert!(self.start_vertex() != self.end_vertex());
        shift_range(r, g.vertex_length(self.start_vertex()))
        //if self.start_vertex() != self.end_vertex() {
        //    shift_range(r, g.node(self.start_vertex().node_id).length)
        //} else {
        //    r
        //}
    }
}

//TODO can be heavily optimized (e.g. no maps, sets, etc)
//TODO support other weights -- currently using max length
//Maybe update to pseudo-code from miniasm paper?
pub struct SbSearchParams {
    pub max_length: usize,
    pub max_diff: usize,
    pub max_count: usize,
}

impl SbSearchParams {
    //all usize values should probably default to max values
    //FIXME provide builder
    pub fn unrestricted() -> SbSearchParams {
        SbSearchParams {
            max_length: usize::MAX,
            max_diff: usize::MAX,
            max_count: usize::MAX,
        }
    }
}

pub fn find_superbubble(g: &Graph, v: Vertex, params: &SbSearchParams) -> Option<Superbubble> {
    find_superbubble_subgraph(g, v, params, None)
}

//TODO handle case when first/last vertex have other outgoing/incoming edges
//last vertex case is almost handled
pub fn find_superbubble_subgraph(g: &Graph, s: Vertex, params: &SbSearchParams,
    consider_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Superbubble> {
    if let Some(f) = consider_vertex_f {
        if !f(s) {
            return None;
        }
    };

    let mut bubble = Superbubble {
        start_vertex: s,
        reached_vertices: HashMap::new(),
        end_vertex: None,
    };

    let outgoing_edge_cnt = |v| {
        match consider_vertex_f {
            None => g.outgoing_edge_cnt(v),
            Some(avail) => g.outgoing_edges(v).iter()
                            .filter(|l| avail(l.end)).count(),
        }
    };

    let incoming_edge_cnt = |v| {
        match consider_vertex_f {
            None => g.incoming_edge_cnt(v),
            Some(avail) => g.incoming_edges(v).iter()
                            .filter(|l| avail(l.start)).count(),
        }
    };

    let outgoing_edges = |v| {
        match consider_vertex_f {
            None => g.outgoing_edges(v),
            Some(avail) => g.outgoing_edges(v).iter().copied()
                            .filter(|l| avail(l.end)).collect(),
        }
    };

    let _incoming_edges = |v| {
        match consider_vertex_f {
            None => g.incoming_edges(v),
            Some(avail) => g.incoming_edges(v).iter().copied()
                            .filter(|l| avail(l.start)).collect(),
        }
    };

    if outgoing_edge_cnt(bubble.start_vertex) < 2
        //same check, but excluding loops
        || outgoing_edges(bubble.start_vertex).iter().filter(|l| l.start != l.end).count() < 2 {
        return None;
    }

    debug!("Adding starting vertex {} to stack", g.v_str(bubble.start_vertex));
    //vertices with all incoming edges considered (can be processed)
    let mut can_be_processed: Vec<Vertex> = vec![bubble.start_vertex];
    bubble.reached_vertices.insert(bubble.start_vertex, (0, 0));

    //reached vertices that can't be processed yet
    let mut not_ready_cnt = 0;
    let mut remaining_incoming: HashMap<Vertex, usize> = HashMap::new();

    while !can_be_processed.is_empty() {
        if bubble.reached_vertices.len() > params.max_count {
            return None;
        }

        let v = can_be_processed.pop().unwrap();
        debug!("Adding vertex {} to the bubble", g.v_str(v));

        if outgoing_edge_cnt(v) == 0 {
            debug!("Hit dead-end");
            return None;
        }

        debug!("Looking at neighbors");
        for l in outgoing_edges(v) {
            let w = l.end;
            if w == bubble.start_vertex {
                return None;
                //FIXME re-enable after dealing with usage wrt start/end symmetry absense
                //if v != self.start_vertex {
                //    //no loops involiving the start vertex
                //    return false;
                //} else {
                //    //unless self-loop
                //    continue;
                //}
            }

            if !bubble.reached_vertices.contains_key(&w) {
                if bubble.reached_vertices.contains_key(&w.rc()) {
                    debug!("Reverse-complement vertex {} was already reached",
                        g.v_str(w.rc()));
                    return None;
                }
                not_ready_cnt += 1;
                remaining_incoming.insert(w, incoming_edge_cnt(w));
                bubble.reached_vertices.insert(w, bubble.link_dist_range(l, g).unwrap());
            }
            let rem_inc = remaining_incoming.get_mut(&w).unwrap();
            *rem_inc -= 1;
            //self.reached_vertices.get(w) =
            bubble.reached_vertices.insert(w,
                merge_range(*bubble.reached_vertices.get(&w).unwrap(), bubble.link_dist_range(l, g).unwrap()));

            if *remaining_incoming.get(&w).unwrap() == 0 {
                can_be_processed.push(w);
                not_ready_cnt -= 1;
            }
        }

        if can_be_processed.len() == 1 && not_ready_cnt == 0 {
            //FIXME second case is not a classic one, check that it works!
            //Also needs more work to get final vertex!!!
            //|| (can_be_processed.len() == 0 && not_ready_cnt == 1)
            //process last vertex?
            let t = can_be_processed.pop().unwrap();
            debug!("End node found! Vertex {}", g.v_str(t));

            let &(min_len, max_len) = bubble.reached_vertices.get(&t).unwrap();

            let v_len = g.vertex_length(t);

            //FIXME it seems like only start_pos is ever checked
            if min_len > v_len && (min_len - v_len) > params.max_length {
                debug!("Length of minimal additional sequence {} exceeded limit {}",
                    min_len - v_len, params.max_length);
                return None;
            }
            if max_len - min_len > params.max_diff {
                debug!("Minimal and maximal lengths differed by {} exceeded limit {}",
                    max_len - min_len, params.max_diff);
                return None;
            }
            bubble.end_vertex = Some(t);
            return Some(bubble);
        }
    }

    debug!("No more nodes could be added");
    debug!("Finished search for starting vertex {}", g.v_str(bubble.start_vertex));
    None
}

pub fn find_all_outer(g: &Graph, params: &SbSearchParams) -> Vec<Superbubble> {
    let mut used_starts = HashSet::new();
    let mut start_2_bubble = HashMap::new();
    for v in g.all_vertices() {
        if used_starts.contains(&v) {
            continue;
        }
        if let Some(bubble) = find_superbubble(g, v, params) {
            //used_starts.insert(bubble.start_vertex());
            used_starts.insert(bubble.end_vertex().rc());
            assert!(!start_2_bubble.contains_key(&bubble.end_vertex().rc()));
            for &w in bubble.inner_vertices() {
                used_starts.insert(w);
                used_starts.insert(w.rc());
                start_2_bubble.remove(&w);
                start_2_bubble.remove(&w.rc());
            }
            start_2_bubble.insert(v, bubble);
        }
    }
    start_2_bubble.into_values().collect()
}

pub type BubbleChain = Vec<Superbubble>;

//TODO maybe switch to Option?
pub fn find_chain_ahead(g: &Graph, init_v: Vertex, params: &SbSearchParams) -> BubbleChain {
    let mut chain = Vec::new();
    //FIXME no need to check here, since we are marking everything, but useful for general code
    let mut v = init_v;

    loop {
        match find_superbubble(g, v, params) {
            None => break,
            Some(bubble) => {
                v = bubble.end_vertex();
                chain.push(bubble);
                if v == init_v {
                    break;
                }
            }
        }
    }
    chain
}

//TODO test
pub fn find_maximal_chain(g: &Graph, mut init_v: Vertex, params: &SbSearchParams) -> BubbleChain {
    let chain_back = find_chain_ahead(g, init_v.rc(), params);
    if !chain_back.is_empty() {
        init_v = chain_back.last().unwrap().end_vertex().rc();
    }
    find_chain_ahead(g, init_v, params)
}

pub fn find_maximal_chains(g: &Graph, params: &SbSearchParams) -> Vec<BubbleChain> {
    let mut considered_start_nodes = HashSet::new();
    let mut maximal_chains = Vec::new();
    for outer_bubble in find_all_outer(g, params) {
        let v = outer_bubble.start_vertex();
        if considered_start_nodes.contains(&v.node_id) {
            continue;
        }
        let chain = find_maximal_chain(g, v, params);
        assert!(!chain.is_empty());
        for bubble in &chain {
            considered_start_nodes.insert(bubble.start_vertex().node_id);
            considered_start_nodes.insert(bubble.end_vertex().node_id);
        }
        maximal_chains.push(chain);
    }
    maximal_chains
}

//will need adjustment if ever 'start' can be same as 'end' in superbubble
pub fn length_range(chain: &[Superbubble], g: &Graph) -> DistRange {
    let mut tot_min = 0;
    let mut tot_max = 0;
    for bubble in chain {
        //TODO implement via negative shift and tuple addition
        let (min, max) = bubble.length_range(g);
        let s_l = g.vertex_length(bubble.start_vertex());
        tot_min += min - s_l;
        tot_max += max - s_l;
    }
    if !chain.is_empty()
        && chain[0].start_vertex() != chain.last().unwrap().end_vertex() {
        let s_l = g.vertex_length(chain[0].start_vertex());
        (tot_min + s_l, tot_max + s_l)
    } else {
        (tot_min, tot_max)
    }
}

//TODO make chain its own structure not to allow empty chains
pub fn longest_path(chain: &[Superbubble], g: &Graph) -> Option<Path> {
    if chain.is_empty() {
        return None;
    }
    let start_vertex = chain[0].start_vertex();
    let mut total = chain[0].longest_path(g);
    for (i, bubble) in chain.iter().enumerate() {
        if i == 0 {
            continue;
        }
        let mut p = bubble.longest_path(g);
        if i == (chain.len() - 1) && bubble.end_vertex() == start_vertex {
            p.trim(1);
        }
        total.extend(p);
    }
    Some(total)
}

pub fn linear_frac(chain: &[Superbubble], g: &Graph) -> f32 {
    assert!(!chain.is_empty());
    let start_vertex = chain[0].start_vertex();
    let mut total_linear = g.vertex_length(start_vertex);
    for (i, bubble) in chain.iter().enumerate() {
        if bubble.end_vertex() != start_vertex {
            total_linear += g.vertex_length(bubble.end_vertex());
        } else {
            assert!(i == chain.len() - 1);
        }
    }
    let longest_path_len = length_range(chain, g).1;
    if total_linear > longest_path_len {
        1.
    } else {
        total_linear as f32 / longest_path_len as f32
    }
}

//FIXME implement flattened vertex iterator even if it has duplicates
pub fn check_chain<F>(chain: &[Superbubble], mut f: F) -> bool
where
F: FnMut(&Vertex) -> bool {
    chain.iter().flat_map(|b| b.vertices()).all(&mut f)
}
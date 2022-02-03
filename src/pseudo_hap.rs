use crate::graph::*;
use crate::graph_algos::*;
use log::{debug, trace};
use std::collections::HashSet;

pub struct LinearBlock {
    instance_path: Path,
    known_alt_nodes: HashSet<usize>,
}

impl LinearBlock {

    pub fn print(&self, g: &Graph) -> String {
        format!("<Block: path={}; known_alts=[{}]>", self.instance_path().print(g),
            self.known_alt_nodes.iter().map(|&node_id| g.name(node_id)).collect::<Vec<&str>>().join(","))
    }

    pub fn instance_path(&self) -> &Path {
        &self.instance_path
    }

    pub fn known_alt_nodes(&self) -> &HashSet<usize> {
        &self.known_alt_nodes
    }

    pub fn all_nodes(&self) -> impl Iterator<Item=usize> + '_ {
        self.instance_path.vertices()
            .iter()
            .map(|v| v.node_id)
            .chain(self.known_alt_nodes.iter().copied())
    }

    fn can_merge_in(&self, other: &LinearBlock) -> bool {
        self.instance_path.can_merge_in(&other.instance_path)
            && other.all_nodes()
             .all(|n| !self.known_alt_nodes.contains(&n))
    }

    fn merge_in(&mut self, other: LinearBlock) {
        debug_assert!(self.can_merge_in(&other));
        self.instance_path.merge_in(other.instance_path);
        self.known_alt_nodes.extend(other.known_alt_nodes.into_iter());
    }

    fn try_merge_in(mut self, other: LinearBlock) -> Option<LinearBlock> {
        if self.can_merge_in(&other) {
            self.merge_in(other);
            Some(self)
        } else {
            None
        }
    }

    fn from_path(path: Path, iter: impl Iterator<Item=Vertex>) -> LinearBlock {
        LinearBlock {
            instance_path: path,
            known_alt_nodes: iter.map(|v| v.node_id).collect(),
        }
    }

    fn from_bubble(g: &Graph, bubble: superbubble::Superbubble) -> LinearBlock {
        let p = bubble.longest_path(g);
        let mut nodes: HashSet<usize> = bubble.vertices().map(|v| v.node_id).collect();
        for v in p.vertices() {
            nodes.remove(&v.node_id);
        }
        LinearBlock {
            instance_path: p,
            known_alt_nodes: nodes.into_iter().collect(),
        }
    }

    fn from_bubble_chain(g: &Graph, bubble_chain: superbubble::BubbleChain) -> LinearBlock {
        assert!(bubble_chain.len() > 0);
        let mut block = Self::vertex_block(bubble_chain[0].start_vertex());
        for b in bubble_chain.into_iter() {
            let b_lb = Self::from_bubble(g, b);
            assert!(block.can_merge_in(&b_lb));
            block.merge_in(b_lb);
        }
        block
    }

    fn vertex_block(v: Vertex) -> LinearBlock {
        LinearBlock {
            instance_path: Path::new(v),
            known_alt_nodes: HashSet::new(),
        }
    }

    fn search_ahead(g: &Graph, v: Vertex, params: &superbubble::SbSearchParams) -> LinearBlock {
        let chain = superbubble::find_chain_ahead(g, v, params);
        if chain.len() > 0 {
            Self::from_bubble_chain(g, chain)
        } else {
            Self::vertex_block(v)
        }
    }

    //fn is_bridge(&self, g: &Graph) -> bool {
    //    g.incoming_edge_cnt(self.instance_path.start()) == 1
    //        && g.outgoing_edge_cnt(self.instance_path.end()) == 1
    //}

    fn reverse_complement(self) -> LinearBlock {
        LinearBlock {
            instance_path: self.instance_path.reverse_complement(),
            known_alt_nodes: self.known_alt_nodes,
            //..self
        }
    }
}

//todo maybe support blocks here? (use block search and is_block method)
fn bridged_by_vertex(g: &Graph, v: Vertex) -> Option<Path> {
    if g.incoming_edge_cnt(v) == 1 && g.outgoing_edge_cnt(v) == 1 {
        let u = g.incoming_edges(v)[0].start;
        let w = g.outgoing_edges(v)[0].end;
        if u.node_id == v.node_id
        || w.node_id == v.node_id
        || w.node_id == u.node_id {
            return None;
        }
        let mut p = Path::from_link(g.incoming_edges(v)[0]);
        p.append(g.outgoing_edges(v)[0]);
        Some(p)
    } else {
        None
    }
}

fn other_outgoing(g: &Graph, v: Vertex, l: Link) -> Option<Link> {
    if g.outgoing_edge_cnt(v) == 2 {
        let alt = g.outgoing_edges(v).iter().copied()
                    .filter(|&x| x != l)
                    .next().unwrap();
        assert!(alt.end != l.end);
        return Some(alt);
    }
    None
}

fn other_incoming(g: &Graph, v: Vertex, l: Link) -> Option<Link> {
    if g.incoming_edge_cnt(v) == 2 {
        let alt = g.incoming_edges(v).iter().copied()
                    .filter(|&x| x != l)
                    .next().unwrap();
        assert!(alt.start != l.start);
        return Some(alt);
    }
    None
}

//Bridge is a path of length 3 with middle vertex having single incoming and single outgoing link
//Returning both middle vertex and entire path
fn bridge_ahead(g: &Graph, v: Vertex) -> Option<Path> {
    let bridges: Vec<Path> = g.outgoing_edges(v).iter()
        .filter_map(|l| bridged_by_vertex(g, l.end))
        .collect();
    if bridges.len() == 1 {
        Some(bridges.into_iter().next().unwrap())
    } else {
        None
    }
}

//TODO move into PrimaryDecomposer and parameterize with superbubble search params
fn unique_block_ahead(g: &Graph, v: Vertex, unique_block_len: usize)
    -> Option<LinearBlock> {
    let block = LinearBlock::search_ahead(g, v, &superbubble::SbSearchParams::unrestricted());
    if block.instance_path.total_length(g) >= unique_block_len {
        Some(block)
    } else {
        None
    }
}

fn unambiguous_outgoing(g: &Graph, v: Vertex) -> Option<Link> {
    match g.outgoing_edge_cnt(v) {
        1 => Some(g.outgoing_edges(v)[0]),
        _ => None,
    }
}

fn forward_extension(g: &Graph, v: Vertex, unique_block_len: usize) -> Option<LinearBlock> {
    //TODO refactor
    extension_via_bridge(g, v, unique_block_len)
        .or(extension_in_deadend(g, v, unique_block_len))
        .or(extension_out_deadend(g, v, unique_block_len))
}

//  x a (for 'alt')
//     \
//- v - w -
fn extension_in_deadend(g: &Graph, v: Vertex, unique_block_len: usize)
-> Option<LinearBlock> {
    let l = unambiguous_outgoing(g, v)?;
    let w = l.end;
    let a = other_incoming(g, w, l)?.start;

    if is_deadend(g, a) {
        let ext_block = LinearBlock::from_path(Path::from_link(l), std::iter::once(a));
        let ext_block = ext_block.try_merge_in(unique_block_ahead(g, w, unique_block_len)?)?;
        Some(ext_block)
    } else {
        None
    }
}

//    a x            a x
//   /       or     /
//- v - w -      - v - o x
//l -- 'horizontal' link
fn extension_out_deadend(g: &Graph, v: Vertex, unique_block_len: usize)
-> Option<LinearBlock> {
    if g.outgoing_edge_cnt(v) == 2 {
        //TODO generalize?
        let mut deadend_links : Vec<Link> = g.outgoing_edges(v).into_iter()
                                        .filter(|&l| is_deadend(g, l.end)).collect();
        deadend_links.sort_by_key(|&l| g.vertex_length(l.end));
        match deadend_links.len() {
            2 => {
                assert!(deadend_links[0].end != deadend_links[1].end);
                let a = deadend_links[0].end;
                let l = deadend_links[1];
                let ext = LinearBlock::from_path(Path::from_link(l), std::iter::once(a));
                return Some(ext);
            }
            1 => {
                let a = deadend_links[0].end;
                let l = other_outgoing(g, v, deadend_links[0]).unwrap();
                let mut ext = LinearBlock::from_path(Path::from_link(l), std::iter::once(a));
                ext.merge_in(unique_block_ahead(g, l.end, unique_block_len)?);
                return Some(ext);
            }
            x => assert!(x == 0),
        }
    }
    None
}

//    s   t
//   /     \
//- u - v - w -
fn extension_via_bridge(g: &Graph, u: Vertex, unique_block_len: usize) -> Option<LinearBlock> {
    if let Some(bridge_p) = bridge_ahead(g, u) {
        assert!(bridge_p.len() == 3);
        //let v = bridge_p.vertices()[1];
        let w = bridge_p.end();
        let s = other_outgoing(g, u, bridge_p.link_at(0))?.end;
        let t = other_incoming(g, w, bridge_p.link_at(1))?.start;

        let ext_block = LinearBlock::from_path(bridge_p,
            admissible_alt_class(g, s, t, unique_block_len)?
            .into_iter());
        let ext_block = ext_block.try_merge_in(unique_block_ahead(g, w, unique_block_len)?)?;
        Some(ext_block)
    } else {
        None
    }
}

//checks if s & t belong to one of considered alt cases and returns alt vertices
fn admissible_alt_class(g: &Graph, s: Vertex, t: Vertex, unique_block_len: usize)
    -> Option<Vec<Vertex>> {
    if s == t {
        //FIXME in this case there can be loop on top of s which won't be added to alt
        return Some(vec![s]);
    }
    if is_deadend(g, s) && is_deadend(g, t) {
        return Some(vec![s, t]);
    }
    joining_vertices(g, s, t, unique_block_len)
}

//returns all vertices lying on
//fn has_alt_path(g: &Graph, s: Vertex, t: Vertex, node_len_thr: usize)
//-> Option<HashSet<Vertex>> {
//    let (visited, _) = graph_algos::bounded_dfs(g, w, node_len_thr);
//    assert!(visited.contains(&w));
//    if visited.contains(&p.end()) {
//        return true;
//    }
//}

//TODO Generalize maybe support simple blocks and/or extra dead-ends (need to then return subgraph info)
fn is_deadend(g: &Graph, v: Vertex) -> bool {
    g.outgoing_edge_cnt(v) == 0 || g.incoming_edge_cnt(v) == 0
}

fn visited_if_reachable(g: &Graph, v: Vertex, w: Vertex,
    direction: dfs::TraversalDirection, max_node_len: usize)
    -> Option<HashSet<Vertex>> {
    let mut dfs = dfs::DFS::new(g, direction);
    dfs.set_max_node_len(max_node_len);
    dfs.extend_blocked(std::iter::once(w));
    dfs.run_from(v);
    if dfs.boundary().contains(&w) {
        Some(dfs.exit_order().iter().copied().collect())
    } else {
        None
    }
}

fn joining_vertices(g: &Graph, s: Vertex, t: Vertex, max_node_len: usize)
    -> Option<Vec<Vertex>> {
    let visited_fwd = visited_if_reachable(g, s, t, dfs::TraversalDirection::FORWARD, max_node_len)?;
    let visited_rev = visited_if_reachable(g, t, s, dfs::TraversalDirection::REVERSE, max_node_len).unwrap();
    let mut reachable: Vec<Vertex> = visited_fwd.intersection(&visited_rev).copied().collect();
    reachable.push(s);
    reachable.push(t);
    Some(reachable)
}

struct PrimaryDecomposer<'a> {
    g: &'a Graph,
    unique_block_len: usize,
    used_nodes: HashSet<usize>,
}

//TODO extend to situations when no single end vertex
//(i.e. blocks ending with simple bubbles)
fn end_vertex(b: &LinearBlock) -> Vertex {
    b.instance_path.end()
}

impl <'a> PrimaryDecomposer<'a> {

    fn new(g: &Graph, unique_block_len: usize) -> PrimaryDecomposer {
        PrimaryDecomposer {
            g,
            unique_block_len,
            used_nodes: HashSet::new(),
        }
    }

    fn extend_forward(&self, block: &mut LinearBlock) -> bool {
        let v = end_vertex(block);
        if let Some(ext) = forward_extension(self.g, v, self.unique_block_len) {
            if ext.all_nodes().all(|n| !self.used_nodes.contains(&n))
               && block.can_merge_in(&ext) {
                block.merge_in(ext);
                return true;
            }
        }
        false
    }

    fn max_extend_forward(&self, block: &mut LinearBlock) -> bool {
        let mut extended = false;
        while self.extend_forward(block) {
            extended = true;
        }
        extended
    }

    //return none if failed to extend
    //FIXME make logic less surprising
    fn extended_block(&self, mut block: LinearBlock) -> Option<LinearBlock> {
        let mut extended = self.max_extend_forward(&mut block);
        let mut rc_block = block.reverse_complement();
        extended |= self.max_extend_forward(&mut rc_block);
        if extended {
            Some(rc_block.reverse_complement())
        } else {
            None
        }
    }

    fn run(&mut self) -> Vec<LinearBlock> {
        let mut resulting_blocks = Vec::new();
        for simple_block in simple_unique_blocks(self.g, self.unique_block_len) {
            if simple_block.all_nodes().all(|n| !self.used_nodes.contains(&n)) {
                if let Some(block) = self.extended_block(simple_block) {
                    assert!(block.all_nodes().all(|n| !self.used_nodes.contains(&n)));
                    self.used_nodes.extend(block.all_nodes());
                    resulting_blocks.push(block);
                }
            }
        }

        for simple_block in simple_unique_blocks(self.g, self.unique_block_len) {
            if simple_block.all_nodes().any(|n| self.used_nodes.contains(&n)) {
                assert!(simple_block.all_nodes().all(|n| self.used_nodes.contains(&n)));
            } else {
                resulting_blocks.push(simple_block);
            }
        }

        resulting_blocks
    }

}

//prioritization step is cheap
fn simple_unique_blocks(g: &Graph, unique_block_len: usize) -> Vec<LinearBlock> {
    use superbubble::*;
    let nodes_in_sccs = scc::nodes_in_sccs(g, &scc::strongly_connected(g));
    let mut used_nodes = HashSet::new();
    //block and it's 'linear fraction' --
    // for single node always 1,
    // for bubble chains -- total fraction of 'joins' in instance paths
    let mut unique_blocks = Vec::new();

    //pub fn linear_frac(chain: &BubbleChain, g: &Graph) -> f32 {
    for chain in find_maximal_chains(g, &SbSearchParams::unrestricted())
                    .into_iter()
                    .filter(|c| check_chain(c, |v| !nodes_in_sccs.contains(&v.node_id))
                                //FIXME think of supporting looped bubble chains
                                && c.first().unwrap().start_vertex() != c.last().unwrap().end_vertex()
                                //FIXME think if makes more sense to check shortest one
                                //but longest path used elsewhere
                                && length_range(c, g).1 >= unique_block_len) {
        assert!(!chain.is_empty());
        assert!(check_chain(&chain, |v| !used_nodes.contains(&v.node_id)));
        for bubble in &chain {
            used_nodes.extend(bubble.vertices().map(|&v| v.node_id));
        }
        let linear_frac = linear_frac(&chain, g);
        unique_blocks.push((LinearBlock::from_bubble_chain(g, chain), linear_frac));
    }

    for (node_id, node) in g.all_nodes().enumerate() {
        if node.length >= unique_block_len && !used_nodes.contains(&node_id) {
            unique_blocks.push((LinearBlock::vertex_block(Vertex::forward(node_id)), 1.));
        }
    }

    unique_blocks.sort_by_cached_key(|(block, lin_frac)| {
        let length = block.instance_path.total_length(g);
        //less linear fraction is better, pulled in one of 11 buckets [0..10] depending on percentage
        let lin_frac_grade = (lin_frac * 10.).round() as u32;
        //within the same 'linear fraction' bucket the longer block the better
        (lin_frac_grade, usize::MAX - length)
    });

    unique_blocks.into_iter().map(|(block, _)| block).collect()
}

pub fn pseudo_hap_decompose(g: &Graph, unique_block_len: usize) -> Vec<LinearBlock> {
    let mut decomposer = PrimaryDecomposer::new(g, unique_block_len);
    decomposer.run()
}

//    s   t
//   /     \
//- u - v - w -
pub fn detect_gap(g: &Graph, u: Vertex) -> Option<GapInfo> {
    if let Some(bridge_p) = bridge_ahead(g, u) {
        assert!(bridge_p.len() == 3);
        //let v = bridge_p.vertices()[1];
        let w = bridge_p.end();
        let s_l = other_outgoing(g, u, bridge_p.link_at(0))?;
        let t_l = other_incoming(g, w, bridge_p.link_at(1))?;
        let s = s_l.end;
        let t = t_l.start;

        if is_deadend(g, s) && is_deadend(g, t) {
            return Some(GapInfo {
                start: s,
                end: t,
                gap_size: (bridge_p.total_length(g) as i64
                           - Path::from_link(s_l).total_length(g) as i64
                           - Path::from_link(t_l).total_length(g) as i64)
            });
        }
    }
    None
}
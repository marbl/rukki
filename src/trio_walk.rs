use crate::graph::*;
use crate::trio::*;
use crate::graph_algos::*;
//FIXME move to common
use scc::only_or_none;
use log::{debug,warn};
use std::collections::{HashSet,HashMap};

// pub fn reachable_ahead(g: &Graph, v: Vertex, node_len_thr: usize) -> HashSet<Vertex> {
//     let (sinks, mut short_ahead) = sinks_ahead(g, v, node_len_thr);
//     short_ahead.extend(&sinks);
//     short_ahead
// }

// pub fn reachable_behind(g: &Graph, v: Vertex, node_len_thr: usize) -> HashSet<Vertex> {
//     let (sources, mut short_behind) = sources_behind(g, v, node_len_thr);
//     short_behind.extend(&sources);
//     short_behind
// }

//TODO optimize
pub fn reachable_between(g: &Graph, v: Vertex, w: Vertex,
    node_len_thr: usize, subgraph_f: Option<&dyn Fn(Vertex)->bool>) -> HashSet<Vertex> {
    let (sinks, mut short_ahead) = dfs::sinks_ahead(g, v, node_len_thr, subgraph_f);
    short_ahead.extend(&sinks);
    let (sources, mut short_behind) = dfs::sources_behind(g, w, node_len_thr, subgraph_f);
    short_behind.extend(&sources);

    short_ahead
        .intersection(&short_behind)
        .copied().collect()
}

const MIN_GAP_SIZE: usize = 1000;
const FILLABLE_BUBBLE_LEN: i64 = 50_000;
const FILLABLE_BUBBLE_DIFF: i64 = 200;

//FIXME use iterators
fn considered_extensions(g: &Graph, v: Vertex,
                consider_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Vec<Link> {
    match consider_vertex_f {
        None => g.outgoing_edges(v),
        Some(avail) => g.outgoing_edges(v).iter().copied()
                        .filter(|l| avail(l.end)).collect(),
    }
}

pub struct ExtensionHelper<'a> {
    g: &'a Graph,
    assignments: &'a AssignmentStorage,
    //FIXME rename to 'consider unassigned' or 'allow unassigned'?
    unassigned_compatible: bool,
}

impl <'a> ExtensionHelper<'a> {

    fn compatible_assignment(&self, node_id: usize, target_group: TrioGroup) -> bool {
        match self.assignments.group(node_id) {
            Some(group) => TrioGroup::compatible(group, target_group),
            None => self.unassigned_compatible,
        }
    }

    fn bearable_assignment(&self, node_id: usize) -> bool {
        match self.assignments.group(node_id) {
            Some(TrioGroup::ISSUE) => false,
            None => self.unassigned_compatible,
            _ => true,
        }
    }

    //FIXME try switching to a &Vertex iterator to simplify calls
    fn only_compatible_of_bearable(&self, v_it: impl Iterator<Item=Vertex> + Clone
        , group: TrioGroup) -> Option<Vertex> {
        //if NOT allowed to use unassigned
        // then check that all the vertices are assigned something
        // (other than ISSUE)
        if v_it.clone().all(|v| self.bearable_assignment(v.node_id)) {
            only_or_none(v_it.filter(|v| self.compatible_assignment(v.node_id, group)))
        } else {
            None
        }
    }

    //FIXME code duplication
    fn only_compatible_of_bearable_link(&self, links: &[Link], group:TrioGroup) -> Option<Link> {
        if links.iter().all(|l| self.bearable_assignment(l.end.node_id)) {
            only_or_none(links.iter().copied()
                .filter(|l| self.compatible_assignment(l.end.node_id, group)))
        } else {
            None
        }
    }

    //maybe move to graph or some GraphAlgoHelper?
    fn group_extension(&self, v: Vertex, group: TrioGroup,
                    consider_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Link> {
        //debug!("Looking at (subset of) outgoing edges for {}", self.g.v_str(v));
        let filtered_outgoing = considered_extensions(self.g, v, consider_vertex_f);

        //If only extension then being unassigned is always ok
        //FIXME Probably obsolete with two-step strategy!
        if filtered_outgoing.len() == 1 {
            let l = filtered_outgoing[0];
            if self.assignments.group(l.end.node_id).map_or(true,
                |g| TrioGroup::compatible(g, group)) {
                debug!("Candidate adjacent extension {} (was unique considered)", self.g.v_str(l.end));
                return Some(l);
            }
        }

        let ext = self.only_compatible_of_bearable_link(&filtered_outgoing, group);
        if let Some(l) = ext {
            debug!("Candidate adjacent extension {}", self.g.v_str(l.end));
        }
        ext
    }

    fn find_compatible_source_sink(&self, v: Vertex, group:TrioGroup, long_node_threshold: usize)
    -> Option<(Vertex, Vertex)> {

        let component = dfs::ShortNodeComponent::search_from(self.g, v, long_node_threshold);

        //check that sources/sinks are clearly separated and that all have assignments
        if !component.simple_boundary() {
            return None;
        }

        let compatible_source = self.only_compatible_of_bearable(component.sources.iter().copied(),
                                        group)?;

        let compatible_sink = self.only_compatible_of_bearable(component.sinks.iter().copied(),
                                        group)?;

        return Some((compatible_source, compatible_sink));
    }

}

//TODO add template parameter
pub struct HaploSearcher<'a> {
    g: &'a Graph,
    assignments: &'a AssignmentStorage,
    extension_helper: ExtensionHelper<'a>,
    solid_len: usize,
    //path intersections by homozygous nodes are always allowed
    allow_intersections: bool,
    //FIXME more reasonable configuration
    //0 -- disabled, 1 -- patch paths, 2 -- fill in small bubbles
    ambig_filling_level: usize,
    used: AssignmentStorage,
    in_sccs: HashSet<usize>,
    small_tangle_index: HashMap<Vertex, scc::LocalizedTangle>,
}

type HaploPath = (Path, usize, TrioGroup);

//FIXME review usage of length threshold!
impl <'a> HaploSearcher<'a> {
    pub fn new(g: &'a Graph, assignments: &'a AssignmentStorage,
        solid_len: usize) -> HaploSearcher<'a> {
        let sccs = scc::strongly_connected(g);
        let mut small_tangle_index = HashMap::new();

        //FIXME parameterize component size separately!
        for small_tangle in scc::find_small_localized(g,
                                                      &sccs,
                                                      solid_len * 3) {
            small_tangle_index.insert(small_tangle.entrance.start, small_tangle);
        }

        HaploSearcher {
            g,
            assignments,
            solid_len,
            allow_intersections: false,
            ambig_filling_level: 1,
            used: AssignmentStorage::new(),
            in_sccs: scc::nodes_in_sccs(g, &sccs),
            extension_helper: ExtensionHelper {
                g,
                assignments,
                unassigned_compatible: false,
            },
            small_tangle_index,
        }
    }

    pub fn new_assigning(g: &'a Graph,
        assignments: &'a AssignmentStorage,
        solid_len: usize) -> HaploSearcher<'a> {
        let mut searcher = Self::new(g, assignments, solid_len);
        searcher.allow_intersections = true;
        searcher.ambig_filling_level = 0;
        searcher.extension_helper.unassigned_compatible = true;
        searcher
    }

    pub fn try_fill_bubbles(&mut self) {
        self.ambig_filling_level = 2;
    }

    //FIXME return from find_all
    pub fn used(&self) -> &AssignmentStorage {
        &self.used
    }

    pub fn take_used(self) -> AssignmentStorage {
        self.used
    }

    //TODO maybe use single length threshold?
    pub fn find_all(&mut self) -> Vec<HaploPath> {
        let mut answer = Vec::new();
        let mut nodes: Vec<(usize, &Node)> = self.g.all_nodes().enumerate().collect();
        nodes.sort_by_key(|(_, n)| n.length);

        for (node_id, node) in nodes.into_iter().rev() {
            if self.used.contains(node_id) || self.in_sccs.contains(&node_id) {
                continue;
            }
            //launch from long, definitely assigned nodes
            if node.length >= self.solid_len && self.assignments.is_definite(node_id) {
                let group = self.assignments.get(node_id).unwrap().group;
                let path = self.haplo_path(Vertex::forward(node_id), group);
                self.used.update_all(path.vertices().iter().map(|v| v.node_id), group);
                self.used.get_mut(path.start().node_id).unwrap().info = String::from("path_boundary");
                self.used.get_mut(path.end().node_id).unwrap().info = String::from("path_boundary");
                answer.push((path, node_id, group));
            }
        }
        answer
    }

    fn haplo_path(&self, v: Vertex, group: TrioGroup) -> Path {
        assert!(self.assignments.group(v.node_id) == Some(group));
        let mut path = Path::new(v);
        self.grow_forward(&mut path, group);
        path = path.reverse_complement();
        self.grow_forward(&mut path, group);
        path.reverse_complement()
    }

    fn aimed_grow_ext(&self, v: Vertex, group: TrioGroup) -> Option<Path> {
        assert!(self.g.vertex_length(v) >= self.solid_len);

        let (u, w) = self.extension_helper.find_compatible_source_sink(v, group, self.solid_len)?;
        assert!(u == v);
        assert!(u.node_id != w.node_id);

        if !self.check_available(w.node_id, group) {
            debug!("Next 'target' vertex {} was unavailable", self.g.v_str(w));
            return None;
        }

        debug!("Found next 'target' vertex {}", self.g.v_str(w));

        let mut reachable_vertices =
            reachable_between(self.g, v, w,
                self.solid_len,
                Some(&|x: Vertex| self.unassigned_or_compatible(x.node_id, group)));

        if reachable_vertices.len() == 0 {
            reachable_vertices =
                reachable_between(self.g, v, w,
                    self.solid_len, None);
        }

        let mut p1 = Path::new(v);
        debug!("Constrained forward extension from {}", self.g.v_str(v));
        self.grow_local(&mut p1, group, w, &|x| reachable_vertices.contains(&x));
        if p1.in_path(w.node_id) {
            //TODO think if actually guaranteed
            assert!(p1.end() == w);
            debug!("Constrained forward search led to complete path");
            return Some(p1);
        }

        let mut p2 = Path::new(w.rc());
        debug!("Constrained backward extension from {}", self.g.v_str(w));
        self.grow_local(&mut p2, group, v.rc(), &|x| reachable_vertices.contains(&x.rc()));
        let p2 = p2.reverse_complement();
        if p2.in_path(v.node_id) {
            //TODO think if actually guaranteed
            assert!(p2.start() == v);
            debug!("Constrained backward search led to complete path");
            return Some(p2);
        }

        //use that multiple copies of the node can't be in path
        if let Some(trim_to) = p1.vertices().iter()
                                        .filter(|x| p2.in_path(x.node_id))
                                        .copied().next() {
            debug!("Paths forward and backward overlapped");
            debug!("Trimming path forward to {}", self.g.v_str(trim_to));
            assert!(p1.trim_to(&trim_to));
            p1.trim(1);
            //FIXME switch to debug_assert
            assert!(p1.vertices().iter().filter(|x| p2.in_path(x.node_id)).next().is_none());
        } else if self.ambig_filling_level > 0 {
            debug!("Will try to patch the gap between forward/backward paths");
            //if paths don't overlap -- try linking
            if let Some(link_p) = self.try_link(p1.end(), p2.start(), group) {
                //TODO simplify? Here we know that p1 and p2 don't have common nodes
                if p1.can_merge_in(&link_p) && link_p.can_merge_in(&p2) {
                    debug!("Patch succesful");
                    p1.merge_in(link_p);
                    p1.merge_in(p2);
                    return Some(p1);
                }
            }
        }

        debug!("Putting ambiguous gap between {} and {}", self.g.v_str(p1.end()), self.g.v_str(p2.start()));
        p1.append_general(GeneralizedLink::AMBIG(GapInfo {
            start: p1.end(),
            end: p2.start(),
            //FIXME use something reasonable
            gap_size: MIN_GAP_SIZE as i64,
        }));
        assert!(p1.can_merge_in(&p2));
        p1.merge_in(p2);
        Some(p1)
    }

    fn grow_forward(&self, path: &mut Path, group: TrioGroup) {
        loop {
            self.aimed_grow(path, group);
            if !self.unguided_grow(path, group) {
                debug!("Stopping extension");
                break;
            }
        }
    }

    //returns true if anything was done and false if couldn't extend
    fn aimed_grow(&self, path: &mut Path, group: TrioGroup) -> bool {
        debug!("Initiating 'guided' extension from {}", self.g.v_str(path.end()));
        let mut something_done = false;
        while let Some(ext) = self.aimed_grow_ext(path.end(), group) {
            assert!(ext.vertices().iter().all(|v| self.check_available(v.node_id, group)));
            debug!("Found extension {}", ext.print(self.g));
            if path.can_merge_in(&ext) {
                debug!("Merging in");
                path.merge_in(ext);
                debug!("Will continue 'guided' extension from {}", self.g.v_str(path.end()));
                something_done = true;
            } else {
                warn!("Couldn't merge in guided extension from {}", self.g.v_str(path.end()));
                break;
            }
        }
        something_done
    }

    //returns true if reached long node and false if ended in issue or couldn't extend anymore
    fn unguided_grow(&self, path: &mut Path, group: TrioGroup) -> bool {
        debug!("Initiating unguided extension from {}", self.g.v_str(path.end()));
        while let Some(ext) = self.unguided_next_or_gap(path.end(), group) {
            if self.check_available_append(path, &ext, group) {
                path.merge_in(ext);
                let v = path.end();
                if self.g.vertex_length(v) >= self.solid_len {
                    debug!("Reached long node {}", self.g.v_str(v));
                    return true;
                }
            } else {
                debug!("Had issue growing beyond {}", self.g.v_str(path.end()));
                return false;
            }
        }
        false
    }

    //TODO maybe consume when grow?
    fn unguided_next_or_gap(&self, v: Vertex, group: TrioGroup) -> Option<Path> {
        self.local_next(v, group, None)
            .or_else(|| self.patch_forward(v, group))
            .or_else(|| self.generalized_patch_forward(v, group))
        //debug!("Was able to patch broken bubble and extend by {grow} nodes");
        //debug!("Was able to patch more general broken haplotype and extend by {grow} nodes");
        //grow += self.merge_in(path, self.find_jump_path_ahead(path.end(), group), group);
        //if grow > 0 {
        //    debug!("Was able to jump (and stitch) ahead by {grow} nodes");
        //    tot_grow += grow;
        //    continue;
        //}
        //grow += self.merge_in(path, self.find_gapped_jump_ahead(path, group), group);
        //if grow > 0 {
        //    debug!("Was able to jump (via ambiguous region) ahead by {grow} nodes");
        //    tot_grow += grow;
        //    continue;
        //}
    }

    //FIXME isn't it obsolete with generalized_patch?
    //FIXME rename
    fn patch_forward(&self, v: Vertex, group: TrioGroup) -> Option<Path> {
        if let Some(gap_info) = self.generalized_gap(v, group, 0) {
            let next_node = gap_info.end.node_id;
            assert!(self.assignments.group(next_node) == Some(group));
            debug!("Identified jump across gap to {}", self.g.v_str(gap_info.end));
            Some(Path::from_general_link(GeneralizedLink::GAP(gap_info)))
        } else {
            None
        }
    }

    fn find_unbroken_alt_candidate(&self, v: Vertex, short_node_threshold: usize) -> Option<(Vertex, i64)> {
        //not necessary, but improves 'symmetry'
        assert!(self.g.vertex_length(v) >= short_node_threshold);

        //dead-end case
        if self.g.outgoing_edge_cnt(v) == 0 {
            let component = dfs::ShortNodeComponent::back_from_long(self.g,
                                            v, short_node_threshold);

            //think of maybe relaxing
            if !component.simple_boundary() {
                return None;
            }

            only_or_none(component.sinks.iter().copied().filter(|&s| s != v))
                .map(|alt| (alt, self.g.vertex_length(alt) as i64 - self.g.vertex_length(v) as i64))
        } else if self.g.outgoing_edge_cnt(v) == 1 {
            //haplotype merge-in case
            let alt = self.g.outgoing_edges(v)[0].end;
            Some((alt, self.g.vertex_length(alt) as i64))
        } else {
            None
        }
    }

    //FIXME add debug prints
    //FIXME rename
    //TODO very asymmetric condition :(
    fn generalized_gap(&self, v: Vertex, group: TrioGroup, short_node_threshold: usize) -> Option<GapInfo> {
        //FIXME might be much easier to augment graph with extra 'gap' links after all!
        if self.g.vertex_length(v) < short_node_threshold {
            return None;
        }
        debug!("Trying to find generalized gap from {}", self.g.v_str(v));
        let (alt, curr_gap_est) = self.find_unbroken_alt_candidate(v, short_node_threshold)?;
        debug!("Was able to find alt node {}", self.g.v_str(alt));

        //FIXME also make it work when alt is short!
        if self.unassigned_or_compatible(alt.node_id, group)
            || self.g.vertex_length(alt) < short_node_threshold {
            debug!("Bad alt");
            return None;
        }

        debug!("Searching for short-node component ahead of {}", self.g.v_str(alt));
        let component = dfs::ShortNodeComponent::ahead_from_long(self.g,
                                        alt, short_node_threshold);

        //think of maybe relaxing
        if !component.simple_boundary()
            || !component.sources.iter().all(|x| self.assignments.is_definite(x.node_id)) {
            debug!("Bad component");
            return None;
        }

        if let Some(&w) = only_or_none(component.sources.iter()
                        .filter(|s| self.assignments.group(s.node_id).unwrap() == group)
                        .filter(|&s| self.g.incoming_edge_cnt(*s) == 0)) {
            debug!("Dead-end case");
            //dead-end case
            return Some(GapInfo {
                start: v,
                end: w,
                gap_size: std::cmp::max(curr_gap_est
                        - self.g.vertex_length(w) as i64, MIN_GAP_SIZE as i64),
            });
        } else if component.sources.len() == 1 {
            //haplotype merge-in case
            assert!(component.sources.iter().next() == Some(&alt));
            if !component.has_deadends
                && component.sinks.iter().all(|x| self.assignments.is_definite(x.node_id)) {
                if let Some(&w) = only_or_none(component.sinks.iter()
                                .filter(|s| self.assignments.group(s.node_id).unwrap() == group)) {
                    debug!("In merge-in case");
                    return Some(GapInfo {
                        start: v,
                        end: w,
                        gap_size: std::cmp::max(curr_gap_est
                                , MIN_GAP_SIZE as i64),
                    });
                }
            }
        }

        debug!("Couldn't find generalized gap");
        None
    }

    //FIXME inline
    //FIXME rename
    fn generalized_patch_forward(&self, v: Vertex, group: TrioGroup) -> Option<Path> {
        //FIXME configure
        if let Some(gap_info) = self.generalized_gap(v, group, 200_000) {
            let next_node = gap_info.end.node_id;
            assert!(self.assignments.group(next_node) == Some(group));
            debug!("Identified jump across 'generalized' gap to {}", self.g.v_str(gap_info.end));
            Some(Path::from_general_link(GeneralizedLink::GAP(gap_info)))
        } else {
            None
        }
    }

    fn long_node(&self, node_id: usize) -> bool {
        self.g.node(node_id).length >= self.solid_len
    }

    fn check_link_vertex(&self, w: Vertex, group: TrioGroup) -> bool {
        let long_node_ahead = |v: Vertex| {
            assert!(self.g.outgoing_edge_cnt(v) == 1);
            self.long_node(self.g.outgoing_edges(v)[0].end.node_id)
        };

        !self.long_node(w.node_id)
            //this check will never allow to patch with unassigned node
            //&& self.assignments.contains(w.node_id)
            && self.check_available(w.node_id, group)
            && self.g.incoming_edge_cnt(w) == 1
            && self.g.outgoing_edge_cnt(w) == 1
            && (long_node_ahead(w)
                || long_node_ahead(w.rc())
                || self.check_assignment(w.node_id, group))
    }

    fn try_link(&self, u: Vertex, w: Vertex, group: TrioGroup) -> Option<Path> {
        debug!("Trying to link {} and {}", self.g.v_str(u), self.g.v_str(w));
        for l in self.g.outgoing_edges(u) {
            if l.end == w {
                debug!("Found direct link {}", self.g.l_str(l));
                return Some(Path::from_link(l))
            }
        }

        let mut outgoing_edges = self.g.outgoing_edges(u);
        outgoing_edges.sort_by(|a, b| self.g.node(b.end.node_id).coverage
                        .partial_cmp(&self.g.node(a.end.node_id).coverage)
                        .unwrap());

        for l in outgoing_edges {
            let v = l.end;
            //TODO think if checks are reasonable //FIXME think if we should check coverage too
            if self.check_link_vertex(v, group) {
                if let Some(l2) = self.g.connector(v, w) {
                    debug!("Was able to link {} and {} via {}",
                        self.g.v_str(u), self.g.v_str(w), self.g.v_str(v));
                    let mut path = Path::from_link(l);
                    path.append(l2);
                    return Some(path);
                }
            }
        }
        None
    }

    fn trivial_bubble_end(&self, u: Vertex, hint_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Vertex> {
        if self.g.outgoing_edge_cnt(u) < 2 {
            return None;
        }
        let mut opt_w = None;
        for v in considered_extensions(self.g, u, hint_f).iter().map(|l1| l1.end) {
            if self.g.incoming_edge_cnt(v) > 1 {
                return None;
            }
            assert!(self.g.incoming_edge_cnt(v) == 1);
            let l2 = only_or_none(considered_extensions(self.g, v, hint_f).into_iter())?;
            if let Some(w) = opt_w {
                if w != l2.end {
                    return None;
                }
            } else {
                opt_w = Some(l2.end);
            }
        }
        assert!(opt_w.is_some());
        opt_w
    }

    //FIXME think if require v assignment
    //TODO sometimes limiting the search here could help
    fn choose_simple_bubble_side(&self, v: Vertex, group: TrioGroup,
        consider_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Path> {
        if self.ambig_filling_level < 2 {
            return None;
        }

        let len_across = |u, v, w| {
            self.g.vertex_length(v) as i64
                - self.g.connector(u, v).unwrap().overlap as i64
                - self.g.connector(v, w).unwrap().overlap as i64
        };

        let end_cov = |l: &Link| self.g.node(l.end.node_id).coverage;
        let w = self.trivial_bubble_end(v, consider_vertex_f)?;
        let filtered_outgoing: Vec<Link> = considered_extensions(self.g, v, consider_vertex_f).into_iter()
                                        .filter(|l| self.unassigned_or_compatible(l.end.node_id, group)).collect();
        if filtered_outgoing.len() > 0
            && filtered_outgoing.iter().all(|&l| self.g.vertex_length(l.end) < self.solid_len)
            && w.node_id != v.node_id {
            let max_len = filtered_outgoing.iter().map(|l| len_across(v, l.end, w)).max().unwrap();
            let min_len = filtered_outgoing.iter().map(|l| len_across(v, l.end, w)).min().unwrap();
            if filtered_outgoing.len() > 1
                && (max_len > FILLABLE_BUBBLE_LEN
                    || (max_len - min_len) > FILLABLE_BUBBLE_DIFF) {
                return None;
            }
            let ext = filtered_outgoing.into_iter()
                                        .max_by(|a, b| end_cov(a).partial_cmp(&end_cov(b)).unwrap())
                                        .unwrap();
            debug!("Candidate simple bubble side {}", self.g.v_str(ext.end));
            Some(Path::from_link(ext))
        } else {
            None
        }
    }

    //TODO improve to actually use group
    //FIXME verify logic of max_length threshold for superbubbles
    fn find_bubble_fill_ahead(&self, v: Vertex, _group: TrioGroup,
        consider_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Path> {
        if self.ambig_filling_level < 2 {
            return None;
        }

        use superbubble::SbSearchParams;
        let sb_params = SbSearchParams {
            max_length: FILLABLE_BUBBLE_LEN as usize,
            max_diff: FILLABLE_BUBBLE_DIFF as usize,
            ..SbSearchParams::unrestricted()
        };
        //TODO think of growing within the bubble if possible (ensyre symmetry)
        let bubble = superbubble::find_superbubble_subgraph(self.g, v, &sb_params,
            consider_vertex_f)?;
        if bubble.inner_vertices().any(|&v| self.g.vertex_length(v) >= self.solid_len) {
            return None;
        }
        let p = bubble.longest_path(self.g);
        debug!("Candidate extension by super-bubble bubble fill {}", p.print(self.g));
        Some(p)
    }

    fn find_bubble_jump_ahead(&self, v: Vertex, _group: TrioGroup,
        consider_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Path> {
        use superbubble::SbSearchParams;
        let sb_params = SbSearchParams {
            //TODO think of relaxing a bit
            max_length: self.solid_len,
            ..SbSearchParams::unrestricted()
        };
        //TODO think of growing within the bubble if possible (ensyre symmetry)
        let bubble = superbubble::find_superbubble_subgraph(self.g, v, &sb_params,
            consider_vertex_f)?;
        if bubble.inner_vertices().any(|&v| self.g.vertex_length(v) >= self.solid_len) {
            return None;
        }
        let w = bubble.end_vertex();
        let gap_est = if bubble.length_range(self.g).0
                        > self.g.vertex_length(v) + self.g.vertex_length(w) + MIN_GAP_SIZE {
            bubble.length_range(self.g).0 - self.g.vertex_length(v) - self.g.vertex_length(w)
        } else {
            MIN_GAP_SIZE
        };
        debug!("Candidate across-bubble jump to {}", self.g.v_str(w));
        Some(Path::from_general_link(GeneralizedLink::AMBIG(GapInfo {
            start: v,
            end: w,
            gap_size: gap_est as i64,
        })))
    }

    fn find_small_tangle_jump_ahead(&self, v: Vertex, _group: TrioGroup) -> Option<Path> {
        let small_tangle = self.small_tangle_index.get(&v)?;
        debug!("Candidate tangle jump to {}", self.g.v_str(small_tangle.exit.end));
        Some(Path::from_general_link(GeneralizedLink::AMBIG(GapInfo {
            start: small_tangle.entrance.start,
            end: small_tangle.exit.end,
            //TODO cache estimated size inside tangle
            gap_size: std::cmp::max(scc::estimate_size_no_mult(small_tangle, self.g),
                                    MIN_GAP_SIZE) as i64,
        })))
    }

    fn unassigned_or_compatible(&self, node_id: usize, group: TrioGroup) -> bool {
        if let Some(assign_group) = self.assignments.group(node_id) {
            if TrioGroup::incompatible(assign_group, group) {
                //if target group is incompatible with initial assignment (incl. ISSUE)
                return false;
            }
        }
        true
    }

    //FIXME maybe stop grow process immediately when this fails
    fn check_available(&self, node_id: usize, target_group: TrioGroup) -> bool {
        if self.in_sccs.contains(&node_id) {
            //is part of some non-trivial SCC
            return false;
        }

        if !self.unassigned_or_compatible(node_id, target_group) {
            return false;
        }

        if !self.allow_intersections {
            if let Some(used_group) = self.used.group(node_id) {
                if TrioGroup::incompatible(used_group, target_group) {
                    //node already used in different haplotype
                    if self.long_node(node_id)
                        && self.assignments.group(node_id) != Some(TrioGroup::HOMOZYGOUS) {
                        //FIXME increase this threshold
                        debug!("Can't reuse long node {} (not initially marked as homozygous) in different haplotype",
                            self.g.name(node_id));
                        return false;
                    }
                } else {
                    //node already used within the same haplotype
                    debug!("Tried to reuse node {} twice within the same haplotype: {:?}",
                            self.g.name(node_id), target_group);
                    return false;
                }
            }
        }
        true
    }

    //returns false if hit some issue (self-intersection, node reuse, etc)
    //TODO optimize
    fn check_available_append(&self, path: &Path, ext: &Path, group: TrioGroup) -> bool {
        path.can_merge_in(ext)
            && ext.links().iter()
                    .all(|l| self.check_available(l.end().node_id, group))
    }

    fn local_next(&self, v: Vertex, group: TrioGroup,
                    constraint_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Path> {
        self.find_small_tangle_jump_ahead(v, group)
            .or_else(|| self.extension_helper.group_extension(v, group, constraint_vertex_f)
                .map(|l| Path::from_link(l)))
            .or_else(|| self.choose_simple_bubble_side(v, group, constraint_vertex_f))
            .or_else(|| self.find_bubble_fill_ahead(v, group, constraint_vertex_f))
            .or_else(|| self.find_bubble_jump_ahead(v, group, constraint_vertex_f))
    }

    ////returns false if ended in issue
    //fn grow(&self, path: &mut Path, group: TrioGroup,
    //                //FIXME typedef?
    //                next_link_f: impl Fn(Vertex)->Option<GeneralizedLink>) -> bool {
    //    let mut tot_grow = 0;
    //    while let Some(l) = next_link_f(path.end()) {
    //        if self.merge_in_if_available(path, l, group) {
    //            tot_grow += 1;
    //        } else {
    //            return false;
    //        }
    //    }
    //    true
    //}

    //returns false if ended in issue
    fn grow_local(&self, path: &mut Path, group: TrioGroup,
                    target_v: Vertex, hint_f: &dyn Fn(Vertex)->bool) -> bool {
        debug!("Locally growing from {}", self.g.v_str(path.end()));
        while let Some(ext) = self.local_next(path.end(), group, Some(hint_f))
                       .or_else(|| self.local_next(path.end(), group, None)) {
            assert!(ext.end() == target_v || !ext.in_path(target_v.node_id));
            if self.check_available_append(path, &ext, group) {
                path.merge_in(ext);
                if path.end() == target_v {
                    break;
                }
            } else {
                debug!("Couldn't append extension {} to the path", ext.print(self.g));
                return false;
            }
        }
        true
    }

    fn check_assignment(&self, node_id: usize, target_group: TrioGroup) -> bool {
        if let Some(assign) = self.assignments.get(node_id) {
            if assign.group == target_group {
                return true;
            }
        }
        false
    }

    ////maybe move to graph or some GraphAlgoHelper?
    //fn unambiguous_extension(&self, v: Vertex) -> Option<Link> {
    //    //TODO simplify?
    //    match self.g.outgoing_edge_cnt(v) {
    //        1 => Some(self.g.outgoing_edges(v)[0]),
    //        _ => None,
    //    }
    //}

}

#[cfg(test)]
mod tests {
    use crate::graph;
    use crate::trio;
    use crate::trio_walk;
    use std::fs;
    use log::info;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn scc_loop_jump() {
        init();

        let graph_fn = "tests/test_graphs/scc_tangle.gfa";
        let assignments_fn = "tests/test_graphs/scc_tangle.ann.csv";
        let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
        let assignments = trio::parse_read_assignments(&g, assignments_fn).unwrap();

        let haplo_searcher = trio_walk::HaploSearcher::new(&g, &assignments, 500_000);
        let path = haplo_searcher.haplo_path(graph::Vertex::forward(g.name2id("utig4-2545")), trio::TrioGroup::PATERNAL);
        assert!(path.len() == 2);
        if let graph::GeneralizedLink::AMBIG(ambig) = path.general_link_at(0) {
            assert!(ambig.gap_size > 900_000 && ambig.gap_size < 1_000_000);
        } else {
            panic!();
        }

        assert_eq!(path.print(&g), String::from("utig4-2545+,AMBIG,utig4-648-"));
    }

    #[test]
    fn gap_jump() {
        init();

        let graph_fn = "tests/test_graphs/test_gap.gfa";
        let assignments_fn = "tests/test_graphs/test_gap.ann.csv";
        let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
        let assignments = trio::parse_read_assignments(&g, assignments_fn).unwrap();

        let haplo_searcher = trio_walk::HaploSearcher::new(&g, &assignments, 500_000);
        for node in ["utig4-1322", "utig4-1320", "utig4-947"] {
            info!("Starting from {}", node);
            println!("Print Starting from {}", node);
            let path = haplo_searcher.haplo_path(graph::Vertex::forward(g.name2id(node)), trio::TrioGroup::MATERNAL);

            assert!(path.len() == 4);
            assert_eq!(path.print(&g), String::from("utig4-947+,utig4-1318-,utig4-1320+,GAP,utig4-1322+"));
            if let graph::GeneralizedLink::GAP(gap) = path.general_link_at(2) {
                assert_eq!(gap.gap_size, 36_423i64);
            } else {
                panic!();
            }
        }
    }
}
use crate::graph::*;
use crate::trio::*;
use crate::graph_algos::*;
use crate::graph_algos::only_or_none;
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
    allow_unassigned: bool,
}

impl <'a> ExtensionHelper<'a> {

    fn compatible_assignment(&self, node_id: usize, target_group: TrioGroup) -> bool {
        match self.assignments.group(node_id) {
            Some(group) => TrioGroup::compatible(group, target_group),
            None => self.allow_unassigned,
        }
    }

    fn bearable_assignment(&self, node_id: usize) -> bool {
        match self.assignments.group(node_id) {
            Some(TrioGroup::ISSUE) => false,
            None => self.allow_unassigned,
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

        //If only extension exists it is always ok if it is unassigned
        //FIXME Probably obsolete with two-step strategy!
        if self.g.outgoing_edge_cnt(v) == 1 {
            let l = self.g.outgoing_edges(v)[0];
            if self.assignments.group(l.end.node_id).map_or(true,
                |g| TrioGroup::compatible(g, group)) {
                debug!("Candidate adjacent extension {} (was unique considered)", self.g.v_str(l.end));
                return Some(l);
            }
        }

        //debug!("Looking at (subset of) outgoing edges for {}", self.g.v_str(v));
        let filtered_outgoing = considered_extensions(self.g, v, consider_vertex_f);
        let ext = self.only_compatible_of_bearable_link(&filtered_outgoing, group);
        if let Some(l) = ext {
            debug!("Candidate adjacent extension {}", self.g.v_str(l.end));
        }
        ext
    }

    fn find_compatible_sink(&self, v: Vertex, group: TrioGroup, solid_len: usize)
    -> Option<Vertex> {
        assert!(self.g.vertex_length(v) >= solid_len);
        let component = dfs::ShortNodeComponent::search_from(self.g, v, solid_len);
        assert!(component.sources.contains(&v));
        debug!("Component -- {}", component.print(self.g));
        debug!("Looking for compatible sink and checking uniqueness");

        //special hairpin case
        if component.sources.len() == 2
            && component.sources.iter()
                .all(|x| self.compatible_assignment(x.node_id, group))
            && component.sinks.iter()
                .map(|&x| x.rc())
                .collect::<HashSet<Vertex>>() == component.sources {
            //let it = component.sinks.iter().copied().filter(|&x| x != v.rc());
            //let a = only_or_none(it);
            debug!("Special hairpin case");
            return only_or_none(component.sinks.iter().copied().filter(|&x| x != v.rc()));
        }

        debug!("Identifying source");
        let s = self.only_compatible_of_bearable(component.sources.iter().copied(),
                                group)?;

        debug!("Identifying sink");
        let t = self.only_compatible_of_bearable(component.sinks.iter().copied(),
                                group)?;

        assert!(s == v);
        if s.node_id == t.node_id {
            debug!("Next 'target' node {} was the same as current one", self.g.v_str(t));
            return None;
        }

        return Some(t);
    }

}

#[derive(Copy, Clone, Debug)]
pub struct HaploSearchSettings {
    //configuring node length thresholds
    pub solid_len: usize,
    pub trusted_len: usize,

    //configuring behavior
    //NB: path intersections by homozygous labeled 'solid' nodes are always allowed
    pub allow_solid_intersections: bool,
    pub allow_unassigned: bool,

    //configuring ambiguous region filling
    //FIXME more reasonable configuration
    //0 -- disabled, 1 -- patch paths, 2 -- fill in small bubbles
    pub ambig_filling_level: usize,
    pub fillable_bubble_len: usize,
    pub fillable_bubble_diff: usize,

    //configuring scaffolding insertion
    pub skippable_tangle_size: usize,
    pub min_gap_size: i64,
    pub default_gap_size: i64,
}

impl Default for HaploSearchSettings {
    fn default() -> Self {
        Self {
            solid_len: 500_000,
            trusted_len: 200_000,
            allow_solid_intersections: false,
            allow_unassigned: false,
            ambig_filling_level: 1,
            fillable_bubble_len: 50_000,
            fillable_bubble_diff: 200,
            skippable_tangle_size: 1_000_000,
            min_gap_size: 1000,
            default_gap_size: 5000,
        }
    }
}

impl HaploSearchSettings {
    pub fn assigning_stage_adjusted(&self) -> HaploSearchSettings {
        HaploSearchSettings {
            allow_solid_intersections: true,
            ambig_filling_level: 0,
            allow_unassigned: true,
            ..*self
        }
    }

    pub fn build_searcher<'a>(&self, g: &'a Graph,
        assignments: &'a AssignmentStorage) -> HaploSearcher<'a> {
            HaploSearcher::new(g, assignments, *self)
    }
}

pub struct HaploSearcher<'a> {
    g: &'a Graph,
    assignments: &'a AssignmentStorage,
    extension_helper: ExtensionHelper<'a>,
    settings: HaploSearchSettings,
    used: AssignmentStorage,
    in_sccs: HashSet<usize>,
    small_tangle_index: HashMap<Vertex, scc::LocalizedTangle>,
}

type HaploPath = (Path, usize, TrioGroup);

impl <'a> HaploSearcher<'a> {
    pub fn new(g: &'a Graph, assignments: &'a AssignmentStorage, settings: HaploSearchSettings) -> HaploSearcher<'a> {
        let sccs = scc::strongly_connected(g);
        let mut small_tangle_index = HashMap::new();

        for small_tangle in scc::find_small_localized(g,
                                                      &sccs,
                                                      settings.skippable_tangle_size) {
            small_tangle_index.insert(small_tangle.entrance.start, small_tangle);
        }

        HaploSearcher {
            g,
            assignments,
            settings,
            used: AssignmentStorage::new(),
            in_sccs: scc::nodes_in_sccs(g, &sccs),
            extension_helper: ExtensionHelper {
                g,
                assignments,
                allow_unassigned: settings.allow_unassigned,
            },
            small_tangle_index,
        }
    }

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
            if node.length >= self.settings.solid_len && self.assignments.is_definite(node_id) {
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
        assert!(self.g.vertex_length(v) >= self.settings.solid_len);

        let w = self.extension_helper.find_compatible_sink(v, group, self.settings.solid_len)?;

        if !self.check_available(w.node_id, group) {
            debug!("Next 'target' node {} was unavailable", self.g.v_str(w));
            return None;
        }

        debug!("Found next 'target' vertex {}", self.g.v_str(w));

        let mut reachable_vertices =
            reachable_between(self.g, v, w,
                self.settings.solid_len,
                Some(&|x: Vertex| self.unassigned_or_compatible(x.node_id, group)));

        if reachable_vertices.len() == 0 {
            reachable_vertices =
                reachable_between(self.g, v, w,
                    self.settings.solid_len, None);
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
            debug_assert!(p1.vertices().iter().filter(|x| p2.in_path(x.node_id)).next().is_none());
        } else if self.settings.ambig_filling_level > 0 {
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
            gap_size: self.settings.default_gap_size,
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
                if self.g.vertex_length(v) >= self.settings.solid_len {
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

    fn unguided_next_or_gap(&self, v: Vertex, group: TrioGroup) -> Option<Path> {
        self.local_next(v, group, None)
            .or_else(|| self.gap_patch(v, group, self.settings.trusted_len))
            //FIXME this one might lead to interesting non-trivial issues
            .or_else(|| self.gap_patch(v, group, 0))
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

    //TODO very asymmetric condition :(
    //TODO might be worthwhile to augment graph with extra 'gap' links in future
    fn generalized_gap_ahead(&self, v: Vertex, group: TrioGroup, short_node_len: usize) -> Option<GapInfo> {
        if self.g.vertex_length(v) < short_node_len {
            return None;
        }
        debug!("Trying to find generalized gap from {}", self.g.v_str(v));
        let (alt, curr_gap_est) = self.find_unbroken_alt_candidate(v, short_node_len)?;
        debug!("Was able to find alt node {}", self.g.v_str(alt));

        //FIXME also make it work when alt is short!
        if self.unassigned_or_compatible(alt.node_id, group)
            || self.g.vertex_length(alt) < short_node_len {
            debug!("Bad alt");
            return None;
        }

        debug!("Searching for short-node component ahead of {}", self.g.v_str(alt));
        let component = dfs::ShortNodeComponent::ahead_from_long(self.g,
                                        alt, short_node_len);

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
                        - self.g.vertex_length(w) as i64, self.settings.min_gap_size),
            });
        } else if component.sources.len() == 1 {
            //haplotype merge-in case
            assert!(component.sources.iter().next() == Some(&alt));
            //FIXME more specific orientation in dead-end check
            if !component.has_deadends
                && component.sinks.iter().all(|x| self.assignments.is_definite(x.node_id)) {
                if let Some(&w) = only_or_none(component.sinks.iter()
                                .filter(|s| self.assignments.group(s.node_id).unwrap() == group)) {
                    debug!("In merge-in case");
                    return Some(GapInfo {
                        start: v,
                        end: w,
                        gap_size: std::cmp::max(curr_gap_est
                                , self.settings.min_gap_size),
                    });
                }
            }
        }

        debug!("Couldn't find generalized gap");
        None
    }

    fn gap_patch(&self, v: Vertex, group: TrioGroup, short_node_len: usize) -> Option<Path> {
        let gap_info = self.generalized_gap_ahead(v, group, short_node_len)?;
        let next_node = gap_info.end.node_id;
        assert!(self.assignments.group(next_node) == Some(group));
        debug!("Identified jump across 'generalized' gap to {}", self.g.v_str(gap_info.end));
        Some(Path::from_general_link(GeneralizedLink::GAP(gap_info)))
    }

    fn long_node(&self, node_id: usize) -> bool {
        self.g.node(node_id).length >= self.settings.solid_len
    }

    fn check_link_vertex(&self, w: Vertex, group: TrioGroup) -> bool {
        let long_node_ahead = |v: Vertex| {
            assert!(self.g.outgoing_edge_cnt(v) == 1);
            let node_id = self.g.outgoing_edges(v)[0].end.node_id;
            self.long_node(node_id)
                && self.assignments.group(node_id) == Some(group)
        };

        !self.long_node(w.node_id)
            //this check will never allow to patch with unassigned node
            //&& self.assignments.contains(w.node_id)
            && self.check_available(w.node_id, group)
            && self.g.incoming_edge_cnt(w) == 1
            && self.g.outgoing_edge_cnt(w) == 1
            && (long_node_ahead(w)
                || long_node_ahead(w.rc())
                || self.assignments.group(w.node_id) == Some(group))
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
            //TODO think if checks are reasonable
            //FIXME think if we should check coverage too
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
        let outgoing = considered_extensions(self.g, u, hint_f);
        if outgoing.len() < 2 {
            return None;
        }
        let mut opt_w = None;
        for v in outgoing.iter().map(|l1| l1.end) {
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

    fn choose_trivial_bubble_side(&self, v: Vertex, group: TrioGroup,
        consider_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Path> {
        if self.settings.ambig_filling_level < 2 {
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
            && filtered_outgoing.iter().all(|&l| self.g.vertex_length(l.end) < self.settings.solid_len)
            && w.node_id != v.node_id {
            let max_len = filtered_outgoing.iter().map(|l| len_across(v, l.end, w)).max().unwrap();
            let min_len = filtered_outgoing.iter().map(|l| len_across(v, l.end, w)).min().unwrap();
            if filtered_outgoing.len() > 1
                && (max_len > self.settings.fillable_bubble_len as i64
                    || (max_len - min_len) > self.settings.fillable_bubble_diff as i64) {
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

    fn connecting_path(&self, v1: Vertex, v2: Vertex, v3: Vertex) -> Path {
        let mut path = Path::from_link(self.g.connector(v1, v2).unwrap());
        path.append(self.g.connector(v2, v3).unwrap());
        path
    }

    //TODO improve to actually use group while choosing a path
    fn find_bubble_fill_ahead(&self, v: Vertex, group: TrioGroup,
        consider_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Path> {

        use superbubble::SbSearchParams;
        let sb_params = SbSearchParams {
            //max_length: self.solid_len,
            ..SbSearchParams::unrestricted()
        };

        //TODO think of growing within the bubble if possible (ensyre symmetry)
        let bubble = superbubble::find_superbubble_subgraph(self.g, v, &sb_params,
            consider_vertex_f)?;
        if bubble.inner_vertices().any(|&v| self.g.vertex_length(v) >= self.settings.solid_len) {
            return None;
        }

        let w = bubble.end_vertex();
        assert!(w.node_id != v.node_id);

        let mut direct_connectors: Vec<Vertex>
            = considered_extensions(self.g, v, consider_vertex_f)
                .into_iter()
                .filter_map(|l1| self.g.connector(l1.end, w))
                .map(|l2| l2.start)
                .filter(|tc_v| self.unassigned_or_compatible(tc_v.node_id, group))
                .collect();

        let cov = |x: &Vertex| self.g.node(x.node_id).coverage;
        direct_connectors.sort_by(|a, b| cov(b)
                        .partial_cmp(&cov(a))
                        .unwrap());

        let length_range = bubble.length_range(self.g);

        if self.settings.ambig_filling_level > 1
                && length_range.1 <= self.settings.fillable_bubble_diff + length_range.0
                && length_range.1 <= self.settings.fillable_bubble_len
                    + self.g.vertex_length(v) + self.g.vertex_length(w) {
            if direct_connectors.len() > 0 {
                let p = self.connecting_path(v, direct_connectors[0], w);
                debug!("Candidate extension by super-bubble fill (direct connector) {}", p.print(self.g));
                return Some(p)
            } else {
                let p = bubble.longest_path(self.g);
                debug!("Candidate extension by super-bubble fill (longest path) {}", p.print(self.g));
                return Some(p)
            }
        }

        //check if any satisfies the filling criteria (once added gap won't be filled in)
        if self.settings.ambig_filling_level > 0 {
            for &direct_conn in &direct_connectors {
                if self.check_link_vertex(direct_conn, group) {
                    let p = self.connecting_path(v, direct_conn, w);
                    debug!("Candidate extension by super-bubble fill (link vertex) {}", p.print(self.g));
                    return Some(p);
                }
            }
        }

        None
    }

    fn find_bubble_jump_ahead(&self, v: Vertex, _group: TrioGroup,
        consider_vertex_f: Option<&dyn Fn(Vertex)->bool>) -> Option<Path> {
        use superbubble::SbSearchParams;
        let sb_params = SbSearchParams {
            //TODO think of relaxing a bit
            max_length: self.settings.solid_len,
            ..SbSearchParams::unrestricted()
        };
        //TODO think of growing within the bubble if possible (ensyre symmetry)
        let bubble = superbubble::find_superbubble_subgraph(self.g, v, &sb_params,
            consider_vertex_f)?;
        if bubble.inner_vertices().any(|&v| self.g.vertex_length(v) >= self.settings.solid_len) {
            return None;
        }
        let w = bubble.end_vertex();
        let gap_est = if bubble.length_range(self.g).0
                        > self.g.vertex_length(v) + self.g.vertex_length(w) + self.settings.min_gap_size as usize {
            (bubble.length_range(self.g).0 - self.g.vertex_length(v) - self.g.vertex_length(w)) as i64
        } else {
            self.settings.min_gap_size
        };
        debug!("Candidate across-bubble jump to {}", self.g.v_str(w));
        Some(Path::from_general_link(GeneralizedLink::AMBIG(GapInfo {
            start: v,
            end: w,
            gap_size: gap_est,
        })))
    }

    fn find_small_tangle_jump_ahead(&self, v: Vertex, _group: TrioGroup) -> Option<Path> {
        let small_tangle = self.small_tangle_index.get(&v)?;
        debug!("Candidate tangle jump to {}", self.g.v_str(small_tangle.exit.end));
        Some(Path::from_general_link(GeneralizedLink::AMBIG(GapInfo {
            start: small_tangle.entrance.start,
            end: small_tangle.exit.end,
            //TODO cache estimated size inside tangle
            gap_size: std::cmp::max(scc::estimate_size_no_mult(small_tangle, self.g) as i64,
                                    self.settings.min_gap_size),
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

        if !self.settings.allow_solid_intersections {
            if let Some(used_group) = self.used.group(node_id) {
                if TrioGroup::incompatible(used_group, target_group) {
                    //node already used in different haplotype
                    if self.long_node(node_id)
                        && self.assignments.group(node_id) != Some(TrioGroup::HOMOZYGOUS) {
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
            .or_else(|| self.choose_trivial_bubble_side(v, group, constraint_vertex_f))
            .or_else(|| self.find_bubble_fill_ahead(v, group, constraint_vertex_f))
            .or_else(|| self.find_bubble_jump_ahead(v, group, constraint_vertex_f))
    }

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
        let assignments = trio::parse_node_assignments(&g, assignments_fn).unwrap();

        let haplo_searcher = trio_walk::HaploSearchSettings::default().build_searcher(&g, &assignments);
        let path = haplo_searcher.haplo_path(graph::Vertex::forward(g.name2id("utig4-2545")), trio::TrioGroup::PATERNAL);
        assert!(path.len() == 2);
        if let graph::GeneralizedLink::AMBIG(ambig) = path.general_link_at(0) {
            assert!(ambig.gap_size > 900_000 && ambig.gap_size < 1_000_000);
        } else {
            panic!();
        }
    }

    #[test]
    fn gap_jump() {
        init();

        let graph_fn = "tests/test_graphs/test_gap.gfa";
        let assignments_fn = "tests/test_graphs/test_gap.ann.csv";
        let g = graph::Graph::read(&fs::read_to_string(graph_fn).unwrap());
        let assignments = trio::parse_node_assignments(&g, assignments_fn).unwrap();

        let haplo_searcher = trio_walk::HaploSearchSettings::default().build_searcher(&g, &assignments);
        for node in ["utig4-1322", "utig4-1320", "utig4-947"] {
            info!("Starting from {}", node);
            println!("Print Starting from {}", node);
            let path = haplo_searcher.haplo_path(graph::Vertex::forward(g.name2id(node)), trio::TrioGroup::MATERNAL);

            assert!(path.len() == 4);
            assert_eq!(path.print(&g), String::from("utig4-947+,utig4-1318-,utig4-1320+,[N36423N],utig4-1322+"));
        }
    }
}
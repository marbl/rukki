use crate::graph::*;
use crate::graph_algos::only_or_none;
use crate::graph_algos::*;
use crate::trio::*;
use itertools::Itertools;
use log::{debug, warn};
use std::collections::{HashMap, HashSet};

//FIXME move to dfs.rs
//TODO optimize
pub fn reachable_between(
    g: &Graph,
    v: Vertex,
    w: Vertex,
    node_len_thr: usize,
    node_f: Option<&dyn Fn(Vertex) -> bool>,
) -> HashSet<Vertex> {
    //let (sinks, mut short_ahead) = dfs::sinks_ahead(g, v, node_len_thr, node_f);
    use dfs::*;
    let mut reachable_between: HashSet<Vertex> = HashSet::new();

    let mut fwd_dfs = DFS::new(g, TraversalDirection::FORWARD, node_f);
    fwd_dfs.set_max_node_len(node_len_thr);
    fwd_dfs.extend_blocked(std::iter::once(w));
    //inner_dfs(g, v, node_len_thr, &mut visited, &mut border);
    fwd_dfs.run_from(v);
    if fwd_dfs.boundary().contains(&w) {
        let fwd_visited = fwd_dfs.visited();
        reachable_between.insert(v);
        reachable_between.insert(w);

        let mut bwd_dfs = DFS::new(g, TraversalDirection::REVERSE, node_f);
        bwd_dfs.set_max_node_len(node_len_thr);
        bwd_dfs.extend_blocked(std::iter::once(v));
        bwd_dfs.run_from(w);
        let bwd_visited = bwd_dfs.visited();
        assert!(bwd_dfs.boundary().contains(&v));

        reachable_between.extend(fwd_visited.intersection(&bwd_visited).copied());
    }
    reachable_between
}

//FIXME use iterators
fn considered_extensions(
    g: &Graph,
    v: Vertex,
    consider_vertex_f: Option<&dyn Fn(Vertex) -> bool>,
) -> Vec<Link> {
    match consider_vertex_f {
        None => g.outgoing_edges(v),
        Some(avail) => g
            .outgoing_edges(v)
            .iter()
            .copied()
            .filter(|l| avail(l.end))
            .collect(),
    }
}

pub struct ExtensionHelper<'a> {
    g: &'a Graph,
    assignments: &'a AssignmentStorage,
    allow_unassigned: bool,
}

impl<'a> ExtensionHelper<'a> {
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

    //FIXME refactor out!
    //FIXME try switching to a &Vertex iterator to simplify calls
    fn only_compatible_of_bearable(
        &self,
        v_it: impl Iterator<Item = Vertex> + Clone,
        group: TrioGroup,
    ) -> Option<Vertex> {
        //if NOT allowed to use unassigned
        // then check that all the vertices are assigned something
        // (other than ISSUE)
        if v_it.clone().all(|v| self.bearable_assignment(v.node_id)) {
            //FIXME remove debug
            debug!(
                "{}",
                v_it.clone()
                    .filter(|v| self.compatible_assignment(v.node_id, group))
                    .map(|x| self.g.v_str(x))
                    .join(", ")
            );
            only_or_none(v_it.filter(|v| self.compatible_assignment(v.node_id, group)))
        } else {
            None
        }
    }

    //FIXME code duplication
    fn only_compatible_of_bearable_link(&self, links: &[Link], group: TrioGroup) -> Option<Link> {
        if links
            .iter()
            .all(|l| self.bearable_assignment(l.end.node_id))
        {
            only_or_none(
                links
                    .iter()
                    .copied()
                    .filter(|l| self.compatible_assignment(l.end.node_id, group)),
            )
        } else {
            None
        }
    }

    //maybe move to graph or some GraphAlgoHelper?
    fn group_extension(
        &self,
        v: Vertex,
        group: TrioGroup,
        consider_vertex_f: Option<&dyn Fn(Vertex) -> bool>,
    ) -> Option<Link> {
        //If only extension exists it is always ok if it is unassigned
        //FIXME Probably obsolete with two-step strategy!
        if self.g.outgoing_edge_cnt(v) == 1 {
            let l = self.g.outgoing_edges(v)[0];
            if self
                .assignments
                .group(l.end.node_id)
                .map_or(true, |g| TrioGroup::compatible(g, group))
            {
                debug!("Candidate unambiguous extension {}", self.g.v_str(l.end));
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

    fn find_assigned_ahead(&self, v: Vertex, group: TrioGroup, solid_len: usize) -> Option<Vertex> {
        let check_unassigned = |x: Vertex| self.assignments.get(x.node_id).is_none();
        let mut dfs = dfs::DFS::new(
            self.g,
            dfs::TraversalDirection::FORWARD,
            Some(&check_unassigned),
        );
        dfs.set_max_node_len(solid_len);
        dfs.run_from(v);

        //could be if solid unassigned node is in the boundary
        if dfs.boundary().iter().any(|&x| check_unassigned(x)) {
            return None;
        }

        only_or_none(
            dfs.boundary()
                .iter()
                .filter(|x| self.compatible_assignment(x.node_id, group))
                .copied(),
        )
    }

    fn find_compatible_sink(
        &self,
        v: Vertex,
        group: TrioGroup,
        solid_len: usize,
    ) -> Option<Vertex> {
        assert!(self.g.vertex_length(v) >= solid_len);
        let component = dfs::ShortNodeComponent::search_from(self.g, v, solid_len);
        assert!(component.sources.contains(&v));
        debug!("Component -- {}", component.print(self.g));
        debug!("Looking for compatible sink and checking uniqueness");

        ////in-haplotype hairpin case
        //if component.sources.len() == 2
        //    && component
        //        .sources
        //        .iter()
        //        .all(|x| self.compatible_assignment(x.node_id, group))
        //    && component
        //        .sinks
        //        .iter()
        //        .map(|&x| x.rc())
        //        .collect::<HashSet<Vertex>>()
        //        == component.sources
        //{
        //    //let it = component.sinks.iter().copied().filter(|&x| x != v.rc());
        //    //let a = only_or_none(it);
        //    debug!("Special hairpin case");
        //    return only_or_none(component.sinks.iter().copied().filter(|&x| x != v.rc()));
        //}

        debug!("Identifying sink");
        //excluding self to handle hairpins
        let t = self.only_compatible_of_bearable(
            component.sinks.iter().copied().filter(|&x| x != v.rc()),
            group,
        )?;

        debug!("Identifying source");
        //excluding reverse-complement to handle hairpins
        let s = self.only_compatible_of_bearable(
            component.sources.iter().copied().filter(|&x| x != t.rc()),
            group,
        )?;

        assert!(s == v);
        if s.node_id == t.node_id {
            debug!(
                "Next 'target' node {} was the same as current one",
                self.g.v_str(t)
            );
            return None;
        }

        Some(t)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct HaploSearchSettings {
    //configuring node length thresholds
    pub solid_len: usize,
    pub trusted_len: usize,

    //configuring behavior
    //Allow reuse of nodes within the same haplotype (otherwise prevented),
    // and use of solid unassigned nodes by different haplotypes
    //NB: path intersections by solid 'homozygous' nodes are always allowed
    //NB2: path self-intersections are always prevented
    pub allow_intersections: bool,
    pub allow_unassigned: bool,

    //fill in small bubbles
    pub fill_bubbles: bool,
    pub max_unique_cov: f64,
    pub fillable_bubble_len: usize,
    pub fillable_bubble_diff: usize,
    pub het_fill_bubble_len: usize,
    pub het_fill_bubble_diff: usize,
    pub good_side_cov_gap: f64,

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
            allow_intersections: false,
            allow_unassigned: false,
            fill_bubbles: true,
            max_unique_cov: f64::MAX,
            fillable_bubble_len: 50_000,
            fillable_bubble_diff: 200,
            het_fill_bubble_len: 50_000,
            het_fill_bubble_diff: 200,
            good_side_cov_gap: 5.,
            skippable_tangle_size: 1_000_000,
            min_gap_size: 1000,
            default_gap_size: 5000,
        }
    }
}

impl HaploSearchSettings {
    pub fn assigning_stage_adjusted(&self) -> HaploSearchSettings {
        HaploSearchSettings {
            allow_intersections: true,
            fill_bubbles: false,
            allow_unassigned: true,
            ..*self
        }
    }
}

pub struct HaploSearcher<'a> {
    g: &'a Graph,
    assignments: &'a AssignmentStorage,
    extension_helper: ExtensionHelper<'a>,
    settings: HaploSearchSettings,
    used: AssignmentStorage,
    small_tangle_index: HashMap<Vertex, scc::LocalizedTangle>,
    raw_cnts: Option<&'a HashMap<usize, TrioInfo>>,
}

pub type HaploPath = (Path, usize, TrioGroup);

impl<'a> HaploSearcher<'a> {
    pub fn new(
        g: &'a Graph,
        assignments: &'a AssignmentStorage,
        settings: HaploSearchSettings,
        raw_cnts: Option<&'a HashMap<usize, TrioInfo>>,
    ) -> HaploSearcher<'a> {
        HaploSearcher {
            g,
            assignments,
            settings,
            used: AssignmentStorage::new(),
            extension_helper: ExtensionHelper {
                g,
                assignments,
                allow_unassigned: settings.allow_unassigned,
            },
            small_tangle_index: HashMap::from_iter(
                scc::find_small_localized(
                    g,
                    &scc::strongly_connected(g),
                    settings.skippable_tangle_size,
                )
                .into_iter()
                .map(|s| (s.entrance.start, s)),
            ),
            raw_cnts,
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
        let mut nodes = self.g.all_nodes().enumerate().collect_vec();
        nodes.sort_by_key(|(_, n)| n.length);

        for (node_id, _node) in nodes.into_iter().rev() {
            //launch from long, definitely assigned nodes
            if !self.used.contains(node_id)
                && self.long_node(node_id)
                && self.assignments.is_definite(node_id)
            {
                let group = self.assignments.get(node_id).unwrap().group;
                let path = self.haplo_path(Vertex::forward(node_id), group);
                self.used
                    .update_all(path.vertices().iter().map(|v| v.node_id), group);
                self.used.get_mut(path.start().node_id).unwrap().info =
                    String::from("path_boundary");
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

    fn solid_aimed_step_ext(&self, v: Vertex, group: TrioGroup) -> Option<Path> {
        assert!(self.long_node(v.node_id));

        let w = self
            .extension_helper
            .find_compatible_sink(v, group, self.settings.solid_len)?;

        //Unifying checks
        //if !self.check_available(w.node_id, group) {
        //    debug!("Next 'target' node {} was unavailable", self.g.v_str(w));
        //    return None;
        //}

        debug!("Found next 'target' vertex {}", self.g.v_str(w));

        self.filling_path_between(v, w, group, true)
    }

    fn assigned_aimed_ext(&self, v: Vertex, group: TrioGroup) -> Option<Path> {
        let w = self
            .extension_helper
            .find_assigned_ahead(v, group, self.settings.solid_len)?;

        debug!("Found next 'assigned' vertex {}", self.g.v_str(w));

        //FIXME do we want to allow gaps here?
        self.filling_path_between(v, w, group, false)
    }

    //FIXME extract to some helper
    fn filling_path_between(
        &self,
        v: Vertex,
        w: Vertex,
        group: TrioGroup,
        allow_gaps: bool,
    ) -> Option<Path> {
        let mut reachable_vertices = reachable_between(
            self.g,
            v,
            w,
            self.settings.solid_len,
            Some(&|x: Vertex| self.unassigned_or_compatible(x.node_id, group)),
        );

        if reachable_vertices.is_empty() {
            reachable_vertices = reachable_between(self.g, v, w, self.settings.solid_len, None);
        }

        let mut p1 = Path::new(v);
        debug!("Constrained forward extension from {} to {}", self.g.v_str(v), self.g.v_str(w));
        self.grow_local(&mut p1, group, w, &|x| reachable_vertices.contains(&x));
        if p1.vertices().contains(&w) {
            assert!(p1.end() == w);
            debug!("Constrained forward search led to complete path");
            return Some(p1);
        }

        let mut p2 = Path::new(w.rc());
        debug!("Constrained backward extension from {} to {}", self.g.v_str(w), self.g.v_str(v));
        self.grow_local(&mut p2, group, v.rc(), &|x| {
            reachable_vertices.contains(&x.rc())
        });
        let p2 = p2.reverse_complement();
        if p2.vertices().contains(&v) {
            assert!(p2.start() == v);
            debug!("Constrained backward search led to complete path");
            return Some(p2);
        }

        if !allow_gaps {
            return None;
        }

        //use that multiple copies of the node can't be in path
        if let Some(trim_to) = p1
            .vertices()
            .iter()
            .filter(|x| p2.in_path(x.node_id))
            .copied()
            .next()
        {
            debug!("Paths forward and backward overlapped");
            debug!("Trimming path forward to {}", self.g.v_str(trim_to));
            assert!(p1.trim_to(&trim_to));
            p1.trim(1);
            debug_assert!(!p1.vertices().iter().any(|x| p2.in_path(x.node_id)));
        }

        debug!(
            "Putting ambiguous gap between {} and {}",
            self.g.v_str(p1.end()),
            self.g.v_str(p2.start())
        );
        p1.append_general(GeneralizedLink::GAP(GapInfo {
            start: p1.end(),
            end: p2.start(),
            //FIXME use something reasonable
            gap_size: self.settings.default_gap_size,
            info: String::from("ambig_path"),
        }));
        assert!(p1.can_merge_in(&p2));
        p1.merge_in(p2);
        Some(p1)
    }

    fn grow_forward(&self, path: &mut Path, group: TrioGroup) {
        loop {
            if self.long_node(path.end().node_id) {
                self.solid_aimed_grow(path, group);
            }
            if !self.unguided_grow_to_solid(path, group) {
                debug!("Stopping extension");
                break;
            }
        }
    }

    //Tries to maximally grow the path forward from a solid node, iteratively trying to guess next solid target
    //returns true if anything was done and false if couldn't extend
    fn solid_aimed_grow(&self, path: &mut Path, group: TrioGroup) {
        debug!(
            "Initiating 'guided' extension from {}",
            self.g.v_str(path.end())
        );
        while let Some(ext) = self.solid_aimed_step_ext(path.end(), group) {
            debug!("Found extension {}", ext.print(self.g));
            if self.check_available_append(path, &ext, group) {
                debug!("Merging in");
                path.merge_in(ext);
                debug!(
                    "Will continue 'guided' extension from {}",
                    self.g.v_str(path.end())
                );
            } else {
                warn!(
                    "Couldn't merge in guided extension from {}",
                    self.g.v_str(path.end())
                );
                break;
            }
        }
    }

    //Tries to maximally grow the path forward until the next solid node (without gessing next solid in advance)
    //returns true if reached solid node and false if ended in issue or couldn't extend anymore
    fn unguided_grow_to_solid(&self, path: &mut Path, group: TrioGroup) -> bool {
        debug!(
            "Initiating unguided extension from {}",
            self.g.v_str(path.end())
        );
        //try to make one step ahead, s.a. gap/tangle/bubble and regular extension
        while let Some(ext) = self.unguided_next_or_gap(path.end(), group) {
            if self.check_available_append(path, &ext, group) {
                path.merge_in(ext);
                let v = path.end();
                if self.long_node(v.node_id) {
                    debug!("Reached solid node {}", self.g.v_str(v));
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
            .or_else(|| self.assigned_aimed_ext(v, group))
            .or_else(|| self.gap_patch(v, group, self.settings.trusted_len))
            //FIXME this one might lead to interesting non-trivial issues
            .or_else(|| self.gap_patch(v, group, 0))
    }

    fn find_unbroken_alt_candidate(
        &self,
        v: Vertex,
        short_node_threshold: usize,
    ) -> Option<(Vertex, i64)> {
        //not necessary, but improves 'symmetry'
        assert!(self.g.vertex_length(v) >= short_node_threshold);

        //dead-end case
        if self.g.outgoing_edge_cnt(v) == 0 {
            let component =
                dfs::ShortNodeComponent::back_from_long(self.g, v, short_node_threshold);

            //think of maybe relaxing
            if !component.simple_boundary() {
                return None;
            }

            only_or_none(component.sinks.iter().copied().filter(|&s| s != v)).map(|alt| {
                (
                    alt,
                    self.g.vertex_length(alt) as i64 - self.g.vertex_length(v) as i64,
                )
            })
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
    fn generalized_gap_ahead(
        &self,
        v: Vertex,
        group: TrioGroup,
        short_node_len: usize,
    ) -> Option<GapInfo> {
        if self.g.vertex_length(v) < short_node_len {
            return None;
        }
        debug!("Trying to find generalized gap from {}", self.g.v_str(v));
        let (alt, curr_gap_est) = self.find_unbroken_alt_candidate(v, short_node_len)?;
        debug!("Was able to find alt node {}", self.g.v_str(alt));

        //FIXME also make it work when alt is short!
        if self.unassigned_or_compatible(alt.node_id, group)
            || self.g.vertex_length(alt) < short_node_len
        {
            debug!("Bad alt");
            return None;
        }

        debug!(
            "Searching for short-node component ahead of {}",
            self.g.v_str(alt)
        );
        let component = dfs::ShortNodeComponent::ahead_from_long(self.g, alt, short_node_len);

        //think of maybe relaxing
        if !component.simple_boundary()
            || !component
                .sources
                .iter()
                .all(|x| self.assignments.is_definite(x.node_id))
        {
            debug!("Bad component");
            return None;
        }

        if let Some(&w) = only_or_none(
            component
                .sources
                .iter()
                .filter(|s| self.assignments.group(s.node_id).unwrap() == group)
                .filter(|&s| self.g.incoming_edge_cnt(*s) == 0),
        ) {
            debug!("Dead-end case");
            //dead-end case
            return Some(GapInfo {
                start: v,
                end: w,
                gap_size: std::cmp::max(
                    curr_gap_est - self.g.vertex_length(w) as i64,
                    self.settings.min_gap_size,
                ),
                info: format!("alt-{}", self.g.name(alt.node_id)),
            });
        } else if component.sources.len() == 1 {
            //haplotype merge-in case
            assert!(component.sources.iter().next() == Some(&alt));
            //FIXME more specific orientation in dead-end check
            if !component.has_deadends
                && component
                    .sinks
                    .iter()
                    .all(|x| self.assignments.is_definite(x.node_id))
            {
                if let Some(&w) = only_or_none(
                    component
                        .sinks
                        .iter()
                        .filter(|s| self.assignments.group(s.node_id).unwrap() == group),
                ) {
                    debug!("In merge-in case");
                    return Some(GapInfo {
                        start: v,
                        end: w,
                        gap_size: std::cmp::max(curr_gap_est, self.settings.min_gap_size),
                        info: format!("alt-{}", self.g.name(alt.node_id)),
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
        debug!(
            "Identified jump across 'generalized' gap to {}",
            self.g.v_str(gap_info.end)
        );
        Some(Path::from_general_link(GeneralizedLink::GAP(gap_info)))
    }

    fn long_node(&self, node_id: usize) -> bool {
        self.g.node(node_id).length >= self.settings.solid_len
    }

    fn homozygous_bubble(&self, v: Vertex, w: Vertex) -> bool {
        self.assignments.group(v.node_id) == Some(TrioGroup::HOMOZYGOUS)
            && self.assignments.group(w.node_id) == Some(TrioGroup::HOMOZYGOUS)
    }

    fn bubble_fill_thresholds(&self, v: Vertex, w: Vertex) -> (usize, usize) {
        if self.homozygous_bubble(v, w) {
            (
                self.settings.het_fill_bubble_len,
                self.settings.het_fill_bubble_diff,
            )
        } else {
            (
                self.settings.fillable_bubble_len,
                self.settings.fillable_bubble_diff,
            )
        }
    }

    fn connecting_path(&self, v1: Vertex, v2: Vertex, v3: Vertex) -> Path {
        let mut path = Path::from_link(self.g.connector(v1, v2).unwrap());
        path.append(self.g.connector(v2, v3).unwrap());
        path
    }

    fn raw_marker_excess(&self, v: &Vertex, group: TrioGroup) -> Option<i64> {
        let raw_cnts = self.raw_cnts?;
        let info = raw_cnts.get(&v.node_id)?;
        let mmp = info.mat as i64 - info.pat as i64;
        match group {
            TrioGroup::MATERNAL => Some(mmp),
            TrioGroup::PATERNAL => Some(-mmp),
            _ => panic!(),
        }
    }

    //TODO improve to actually use group while choosing a path
    fn find_bubble_fill_ahead(
        &self,
        v: Vertex,
        group: TrioGroup,
        consider_vertex_f: Option<&dyn Fn(Vertex) -> bool>,
    ) -> Option<Path> {
        //TODO think of growing within the bubble if possible (ensure symmetry)
        let bubble = superbubble::find_superbubble_subgraph(
            self.g,
            v,
            &superbubble::SbSearchParams::unrestricted(),
            consider_vertex_f,
        )?;
        if bubble.inner_vertices().any(|&x| self.long_node(x.node_id)) {
            return None;
        }

        let w = bubble.end_vertex();
        assert!(w.node_id != v.node_id);

        let length_range = bubble.length_range(self.g);

        let (fillable_bubble_len, fillable_bubble_diff) = self.bubble_fill_thresholds(v, w);

        if self.settings.fill_bubbles
            && length_range.1 <= fillable_bubble_diff + length_range.0
            && length_range.1
                <= fillable_bubble_len + self.g.vertex_length(v) + self.g.vertex_length(w)
            && self.bubble_filling_cov_check(v)
            && self.bubble_filling_cov_check(w)
        {
            let cov = |x: &Vertex| self.g.node(x.node_id).coverage;

            //Filling the bubble
            let mut direct_connectors = considered_extensions(self.g, v, consider_vertex_f)
                .into_iter()
                .filter_map(|l1| self.g.connector(l1.end, w))
                .map(|l2| l2.start)
                .filter(|tc_v| self.unassigned_or_compatible(tc_v.node_id, group))
                .collect_vec();

            if !direct_connectors.is_empty() {
                let c = if self.homozygous_bubble(v, w) {
                    let filtered_connectors = direct_connectors
                        .iter()
                        .filter(|&c| !self.used.contains(c.node_id))
                        .filter(|&c| {
                            self.settings.good_side_cov_gap == 0.
                                || self.settings.good_side_cov_gap * cov(c)
                                    > (cov(&v) + cov(&w)) / 2. - 1e-5
                        })
                        .cloned()
                        .collect_vec();

                    if !filtered_connectors.is_empty() {
                        direct_connectors = filtered_connectors;
                    }

                    //direct_connectors.sort_by_key(|c| cov(c));
                    //direct_connectors.into_iter().max_by_key(|&c| (cov(c), self.raw_marker_excess(c, group)));
                    direct_connectors
                        .into_iter()
                        .max_by(|a, b| {
                            (self.raw_marker_excess(a, group).unwrap_or_default(), cov(a))
                                .partial_cmp(&(
                                    self.raw_marker_excess(b, group).unwrap_or_default(),
                                    cov(b),
                                ))
                                .unwrap()
                        })
                        .unwrap()
                } else {
                    direct_connectors
                        .into_iter()
                        .max_by(|a, b| cov(a).partial_cmp(&cov(b)).unwrap())
                        .unwrap()
                };

                let p = self.connecting_path(v, c, w);
                debug!(
                    "Candidate extension by super-bubble fill (direct connector) {}",
                    p.print(self.g)
                );
                Some(p)
            } else {
                let p = bubble.longest_path(self.g);
                debug!(
                    "Candidate extension by super-bubble fill (longest path) {}",
                    p.print(self.g)
                );
                Some(p)
            }
        } else {
            //Jumping across bubble
            let gap_est = if bubble.length_range(self.g).0
                > self.g.vertex_length(v)
                    + self.g.vertex_length(w)
                    + self.settings.min_gap_size as usize
            {
                (bubble.length_range(self.g).0 - self.g.vertex_length(v) - self.g.vertex_length(w))
                    as i64
            } else {
                self.settings.min_gap_size
            };
            debug!("Candidate across-bubble jump to {}", self.g.v_str(w));
            Some(Path::from_general_link(GeneralizedLink::GAP(GapInfo {
                start: v,
                end: w,
                gap_size: gap_est,
                info: String::from("ambig_bubble"),
            })))
        }
    }

    fn find_small_tangle_jump_ahead(&self, v: Vertex, _group: TrioGroup) -> Option<Path> {
        let small_tangle = self.small_tangle_index.get(&v)?;
        debug!(
            "Candidate tangle jump to {}",
            self.g.v_str(small_tangle.exit.end)
        );
        Some(Path::from_general_link(GeneralizedLink::GAP(GapInfo {
            start: small_tangle.entrance.start,
            end: small_tangle.exit.end,
            //TODO cache estimated size inside tangle
            gap_size: std::cmp::max(
                scc::estimate_size_no_mult(small_tangle, self.g) as i64,
                self.settings.min_gap_size,
            ),
            info: String::from("tangle"),
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
        if !self.unassigned_or_compatible(node_id, target_group) {
            return false;
        }

        if !self.settings.allow_intersections {
            if let Some(used_group) = self.used.group(node_id) {
                if TrioGroup::incompatible(used_group, target_group) {
                    //node already used in different haplotype
                    if self.long_node(node_id)
                        && self.assignments.group(node_id) != Some(TrioGroup::HOMOZYGOUS)
                    {
                        assert!(self.assignments.group(node_id).is_none());
                        warn!("Can't reuse long node {} (not initially marked as homozygous) in different haplotype",
                            self.g.name(node_id));
                        return false;
                    }
                } else {
                    //node already used within the same haplotype
                    //TODO consider allowing when deduplication is implemented
                    debug!(
                        "Tried to reuse node {} twice within the same haplotype: {:?}",
                        self.g.name(node_id),
                        target_group
                    );
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
            && ext
                .links()
                .iter()
                .all(|l| self.check_available(l.end().node_id, group))
    }

    fn bubble_filling_cov_check(&self, v: Vertex) -> bool {
        assert!(self.settings.fill_bubbles && self.settings.max_unique_cov >= 0.);
        (self.settings.max_unique_cov > 0.
            && (self.g.node(v.node_id).coverage - 1e-5) < self.settings.max_unique_cov)
            || self.long_node(v.node_id)
            || self.assignments.group(v.node_id) == Some(TrioGroup::HOMOZYGOUS)
    }

    fn local_next(
        &self,
        v: Vertex,
        group: TrioGroup,
        constraint_vertex_f: Option<&dyn Fn(Vertex) -> bool>,
    ) -> Option<Path> {
        self.find_small_tangle_jump_ahead(v, group)
            .or_else(|| {
                self.extension_helper
                    .group_extension(v, group, constraint_vertex_f)
                    .map(Path::from_link)
            })
            .or_else(|| self.find_bubble_fill_ahead(v, group, constraint_vertex_f))
    }

    //returns false if ended in issue
    //FIXME remove return value
    fn grow_local(
        &self,
        path: &mut Path,
        group: TrioGroup,
        target_v: Vertex,
        hint_f: &dyn Fn(Vertex) -> bool,
    ) -> bool {
        debug!("Locally growing from {}", self.g.v_str(path.end()));
        while let Some(ext) = self
            .local_next(path.end(), group, Some(hint_f))
            //FIXME review logic
            .or_else(|| self.local_next(path.end(), group, None))
        {
            assert!(ext.end() == target_v || !ext.in_path(target_v.node_id));
            if self.check_available_append(path, &ext, group) {
                path.merge_in(ext);
                if path.end() == target_v {
                    break;
                }
            } else {
                debug!(
                    "Couldn't append extension {} to the path",
                    ext.print(self.g)
                );
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
    use crate::trio_walk::HaploSearcher;
    use log::info;
    use std::fs;

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

        let haplo_searcher = HaploSearcher::new(
            &g,
            &assignments,
            trio_walk::HaploSearchSettings::default(),
            None,
        );
        let path = haplo_searcher.haplo_path(
            graph::Vertex::forward(g.name2id("utig4-2545")),
            trio::TrioGroup::PATERNAL,
        );
        assert!(path.len() == 2);
        if let graph::GeneralizedLink::GAP(gap) = path.general_link_at(0) {
            assert!(gap.gap_size > 900_000 && gap.gap_size < 1_000_000);
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

        let haplo_searcher = HaploSearcher::new(
            &g,
            &assignments,
            trio_walk::HaploSearchSettings::default(),
            None,
        );
        for node in ["utig4-1322", "utig4-1320", "utig4-947"] {
            info!("Starting from {}", node);
            println!("Print Starting from {node}");
            let path = haplo_searcher.haplo_path(
                graph::Vertex::forward(g.name2id(node)),
                trio::TrioGroup::MATERNAL,
            );

            assert!(path.len() == 4);
            assert_eq!(
                path.print(&g),
                String::from(
                    "utig4-947+,utig4-1318-,utig4-1320+,[N36423N:alt-utig4-1319],utig4-1322+"
                )
            );
        }
    }
}

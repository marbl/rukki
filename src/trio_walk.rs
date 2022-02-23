use crate::graph::*;
use crate::trio::*;
use crate::pseudo_hap::*;
use crate::graph_algos::*;
//FIXME move to common
use scc::only_or_none;
use log::debug;
use std::collections::{HashSet,HashMap};

const MIN_GAP_SIZE: usize = 1000;

pub struct ExtensionHelper<'a> {
    g: &'a Graph,
    assignments: &'a AssignmentStorage<'a>,
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

    //FIXME use iterators
    fn available_extensions(&self, v: Vertex,
                    available_vertex_f: Option<&dyn Fn(&Vertex)->bool>) -> Vec<Link> {
        match available_vertex_f {
            None => self.g.outgoing_edges(v),
            Some(avail) => self.g.outgoing_edges(v).iter().copied()
                            .filter(|l| avail(&l.end)).collect(),
        }
    }

    //maybe move to graph or some GraphAlgoHelper?
    fn group_extension(&self, v: Vertex, group: TrioGroup,
                    available_vertex_f: Option<&dyn Fn(&Vertex)->bool>) -> Option<Link> {
        let filtered_outgoing = self.available_extensions(v, available_vertex_f);

        //If only extension then being unassigned is Ok
        //FIXME Probably obsolete with two-step strategy!
        if filtered_outgoing.len() == 1 {
            let l = filtered_outgoing[0];
            if self.assignments.group(l.end.node_id).map_or(true,
                |g| TrioGroup::compatible(g, group)) {
                return Some(l);
            }
        }

        self.only_compatible_of_bearable_link(&filtered_outgoing, group)
    }

    fn potential_jump_ext(&self, v: Vertex, group: TrioGroup,
        long_node_threshold: usize) -> Option<Vertex> {
        //Currently behavior is quite conservative:
        //1. all long nodes ahead should have assignment
        //2. only one should have correct assignment
        //3. this one should have unambiguous path backward to the vertex maybe stopping one link away
        let (long_ahead, _) = dfs::sinks_ahead(self.g, v, long_node_threshold);

        //long_ahead.retain(|x| x != &v);

        //println!("Long ahead: {}", long_ahead.iter().map(|x| self.g.v_str(*x)).collect::<Vec<String>>().join(";"));

        self.only_compatible_of_bearable(long_ahead.iter().filter(|&x| x != &v).copied(), group)
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
    assignments: &'a AssignmentStorage<'a>,
    extension_helper: ExtensionHelper<'a>,
    long_node_threshold: usize,
    //path intersections by homozygous nodes are always allowed
    allow_intersections: bool,
    used: AssignmentStorage<'a>,
    in_sccs: HashSet<usize>,
    small_tangle_index: HashMap<Vertex, scc::LocalizedTangle>,
}

//FIXME review usage of length threshold!
impl <'a> HaploSearcher<'a> {

    pub fn new(g: &'a Graph, assignments: &'a AssignmentStorage<'a>,
        long_node_threshold: usize) -> HaploSearcher<'a> {
        let sccs = scc::strongly_connected(g);
        let mut small_tangle_index = HashMap::new();

        //FIXME parameterize component size separately!
        for small_tangle in scc::find_small_localized(g,
                                                      &sccs,
                                                      long_node_threshold * 3) {
            small_tangle_index.insert(small_tangle.entrance.start, small_tangle);
        }

        HaploSearcher{
            g,
            assignments,
            long_node_threshold,
            allow_intersections: false,
            used: AssignmentStorage::new(g),
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
        assignments: &'a AssignmentStorage<'a>,
        long_node_threshold: usize) -> HaploSearcher<'a> {
        let mut searcher = Self::new(g, assignments, long_node_threshold);
        searcher.allow_intersections = true;
        searcher.extension_helper.unassigned_compatible = true;
        searcher
    }

    pub fn used(&self) -> &AssignmentStorage<'a> {
        &self.used
    }

    pub fn take_used(self) -> AssignmentStorage<'a> {
        self.used
    }

    //TODO maybe use single length threshold?
    pub fn find_all(&mut self) -> Vec<(Path, usize, TrioGroup)> {
        let mut answer = Vec::new();
        let mut nodes: Vec<(usize, &Node)> = self.g.all_nodes().enumerate().collect();
        nodes.sort_by_key(|(_, n)| n.length);

        for (node_id, node) in nodes.into_iter().rev() {
            if self.used.contains(node_id) || self.in_sccs.contains(&node_id) {
                continue;
            }
            //launch from long, definitely assigned nodes
            if node.length >= self.long_node_threshold && self.assignments.is_definite(node_id) {
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
        self.grow_jump_forward(&mut path, group);
        path = path.reverse_complement();
        self.grow_jump_forward(&mut path, group);
        path.reverse_complement()
    }

    //TODO maybe consume when grow?
    fn grow_jump_forward(&self, path: &mut Path, group: TrioGroup) -> usize {
        let mut tot_grow = 0;
        debug!("Trying to extend forward from vertex {}", self.g.v_str(path.end()));
        loop {
            //FIXME refactor, very ugly!
            let mut grow = self.local_grow(path, group, true, None);
            if grow > 0 {
                debug!("Was able to locally extend by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.jump_forward(path, self.find_jump_path_ahead(path.end(), group), group, None);
            if grow > 0 {
                debug!("Was able to jump (and stitch) ahead by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.jump_forward(path, self.find_gapped_jump_ahead(path, group), group, None);
            if grow > 0 {
                debug!("Was able to jump (via ambiguous region) ahead by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.patch_forward(path, group);
            if grow > 0 {
                debug!("Was able to patch broken bubble and extend by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.generalized_patch_forward(path, group);
            if grow > 0 {
                debug!("Was able to patch more general broken haplotype and extend by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            break;
        }
        tot_grow
    }

    fn patch_forward(&self, path: &mut Path, group: TrioGroup) -> usize {
        let v = path.end();
        if self.g.outgoing_edge_cnt(v) == 0 /* v is dead-end */
            && self.g.incoming_edge_cnt(v) == 1 {
            //FIXME maybe check that the 'joining' node is from other haplotype or is unassigned?
            //todo maybe support case when dead-ends are themselves unassigned? (trivial procedure stops being 'symmetric'))
            if let Some(gap_info) = detect_gap(self.g, self.g.incoming_edges(v)[0].start) {
                let next_node = gap_info.end.node_id;
                if self.assignments.group(next_node) == Some(group)
                    && !path.in_path(next_node)
                    && self.check_available(next_node, group) {
                    path.append_general(GeneralizedLink::GAP(gap_info));
                    return 1;
                }
            }
        }
        0
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
    //TODO very asymmetric condition :(
    fn generalized_gap(&self, v: Vertex, group: TrioGroup, short_node_threshold: usize) -> Option<GapInfo> {
        //FIXME might be much easier to augment graph with extra 'gap' links after all!
        if self.g.vertex_length(v) < short_node_threshold {
            return None;
        }
        let (alt, curr_gap_est) = self.find_unbroken_alt_candidate(v, short_node_threshold)?;
        if self.assignments.is_definite(alt.node_id)
            && TrioGroup::incompatible(group, self.assignments.group(alt.node_id).unwrap()) {
            //debug!("Searching for component ahead from {}", self.g.v_str(alt));
            //FIXME also make it work when alt is short!
            if self.g.vertex_length(alt) < short_node_threshold {
                return None;
            }
            let component = dfs::ShortNodeComponent::ahead_from_long(self.g,
                                            alt, short_node_threshold);

            //think of maybe relaxing
            if !component.simple_boundary() {
                return None;
            }

            if component.sources.iter().all(|x| self.assignments.is_definite(x.node_id)) {
                if let Some(&w) = only_or_none(component.sources.iter()
                                .filter(|s| self.assignments.group(s.node_id).unwrap() == group)
                                .filter(|&s| self.g.incoming_edge_cnt(*s) == 0)) {
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
                            //FIXME code duplication!
                            //debug!("Haplotype merge-in case success");
                            return Some(GapInfo {
                                start: v,
                                end: w,
                                gap_size: std::cmp::max(curr_gap_est
                                        , MIN_GAP_SIZE as i64),
                            });
                        }
                    }
                }
            }
        }
        None
    }

    fn generalized_patch_forward(&self, path: &mut Path, group: TrioGroup) -> usize {
        let v = path.end();
        //FIXME configure
        if let Some(gap_info) = self.generalized_gap(v, group, 100_000) {
            let next_node = gap_info.end.node_id;
            assert!(self.assignments.group(next_node) == Some(group));
            if !path.in_path(next_node)
                && self.check_available(next_node, group) {
                path.append_general(GeneralizedLink::GAP(gap_info));
                return 1;
            }
        }
        0
    }

    fn jump_forward(&self, path: &mut Path, opt_jump: Option<Path>, group: TrioGroup,
        available_vertex_f: Option<&dyn Fn(&Vertex)->bool>) -> usize {
        match opt_jump {
            None => 0,
            Some(jump) => {
                assert!(jump.len() > 1);
                assert!(path.end() == jump.start());
                //FIXME improve logging!
                if path.can_merge_in(&jump)
                    && (&(jump.vertices())[0..(jump.len() - 1)]).iter().all(|v| !self.in_sccs.contains(&v.node_id))
                    && jump.vertices().iter().all(|v| self.check_available(v.node_id, group))
                    && (available_vertex_f.is_none() || jump.vertices().iter().all(|v| available_vertex_f.unwrap()(v))) {
                    let add_on = jump.len() - 1;
                    path.merge_in(jump);
                    add_on
                } else {
                    0
                }
            }
        }
    }

    fn try_link(&self, mut path: Path, v: Vertex) -> Path {
        for l in self.g.outgoing_edges(path.end()) {
            if l.end == v {
                path.append(l);
                break;
            }
        }
        path
    }

    fn long_node(&self, node_id: usize) -> bool {
        self.g.node(node_id).length >= self.long_node_threshold
    }

    fn link_vertex_check(&self, w: Vertex, group: TrioGroup) -> bool {
        let long_node_ahead = |v: Vertex| {
            assert!(self.g.outgoing_edge_cnt(v) == 1);
            self.long_node(self.g.outgoing_edges(v)[0].end.node_id)
        };

        !self.long_node(w.node_id)
            //this check will never allow to patch with unassigned node
            && self.assignments.contains(w.node_id)
            && self.extension_helper.compatible_assignment(w.node_id, group)
            && self.g.incoming_edge_cnt(w) == 1
            && self.g.outgoing_edge_cnt(w) == 1
            && (long_node_ahead(w)
                || long_node_ahead(w.rc())
                || self.check_assignment(w.node_id, group))
    }

    fn try_link_with_vertex(&self, mut path: Path, v: Vertex, group: TrioGroup) -> Path {
        let mut outgoing_edges = self.g.outgoing_edges(path.end());
        outgoing_edges.sort_by(|a, b| self.g.node(b.end.node_id).coverage
                        .partial_cmp(&self.g.node(a.end.node_id).coverage)
                        .unwrap());

        for l in outgoing_edges {
            let w = l.end;
            //TODO think if checks are reasonable //FIXME think if we should check coverage too
            if !path.in_path(w.node_id) && self.link_vertex_check(w, group) {
                if let Some(l2) = self.g.connector(w, v) {
                    debug!("Was able to link {} via {}", self.g.v_str(v), self.g.v_str(w));
                    path.append(l);
                    path.append(l2);
                    break;
                }
            }
        }
        path
    }

    //FIXME inline
    fn reachable_short(&self, v: Vertex, node_len_thr: usize) -> HashSet<Vertex> {
        dfs::sinks_ahead(self.g, v, node_len_thr).1
    }

    fn find_jump_path_ahead(&self, v: Vertex, group: TrioGroup) -> Option<Path> {
        let potential_ext = self.extension_helper.potential_jump_ext(v, group, self.long_node_threshold)?;
        debug!("Unique potential extension {}", self.g.v_str(potential_ext));
        let mut p = Path::new(potential_ext.rc());
        debug!("Growing path forward from {}", self.g.v_str(potential_ext.rc()));
        self.local_grow(&mut p, group, false,
            Some(&|x: &Vertex| {self.reachable_short(v, self.long_node_threshold).contains(&x.rc())}));
        debug!("Found path {}", p.print(self.g));
        if !p.in_path(v.node_id) {
            debug!("Tried linking via vertex");
            p = self.try_link_with_vertex(p, v.rc(), group);
        }
        if !p.in_path(v.node_id) {
            debug!("Tried linking");
            p = self.try_link(p, v.rc());
        }
        if p.trim_to(&v.rc()) {
            assert!(p.len() > 1);
            let p = p.reverse_complement();
            debug!("Successfully found jump, path {}", p.print(self.g));
            return Some(p);
        }
        debug!("Couldn't trim to vertex {}", self.g.v_str(v.rc()));
        None
    }

    fn find_bubble_jump_ahead(&self, v: Vertex, _group: TrioGroup) -> Option<Path> {
        use superbubble::SbSearchParams;
        let sb_params = SbSearchParams {
            //TODO think of relaxing a bit
            max_length: self.long_node_threshold,
            ..SbSearchParams::unrestricted()
        };
        //TODO think of growing within the bubble if possible (ensyre symmetry)
        let bubble = superbubble::find_superbubble(self.g, v, &sb_params)?;
        let w = bubble.end_vertex();
        let gap_est = if bubble.length_range(self.g).0
                        > self.g.vertex_length(v) + self.g.vertex_length(w) + MIN_GAP_SIZE {
            bubble.length_range(self.g).0 - self.g.vertex_length(v) - self.g.vertex_length(w)
        } else {
            MIN_GAP_SIZE
        };
        Some(Path::from_general_link(GeneralizedLink::AMBIG(GapInfo {
            start: v,
            end: w,
            gap_size: gap_est as i64,
        })))
    }

    fn find_small_tangle_jump_ahead(&self, v: Vertex, _group: TrioGroup) -> Option<Path> {
        let small_tangle = self.small_tangle_index.get(&v)?;
        Some(Path::from_general_link(GeneralizedLink::AMBIG(GapInfo {
            start: small_tangle.entrance.start,
            end: small_tangle.exit.end,
            //TODO cache estimated size inside tangle
            gap_size: std::cmp::max(scc::estimate_size_no_mult(small_tangle, self.g),
                                    MIN_GAP_SIZE) as i64,
        })))
    }

    //FIXME rename?
    fn find_gapped_jump_ahead(&self, path: &Path, group:TrioGroup) -> Option<Path> {
        let v = path.end();
        let (u, w) = self.extension_helper.find_compatible_source_sink(v, group, self.long_node_threshold)?;
        if !path.vertices().contains(&u) {
            return None;
        }
        //FIXME optimize, this info should be logged within the short node component
        if !dfs::sinks_ahead(self.g, v, self.long_node_threshold).0.contains(&w) {
            //w can't be reached from v
            return None;
        }
        debug!("Unique potential extension {}", self.g.v_str(w));
        let mut p = Path::new(w.rc());
        debug!("Growing path forward from {}", self.g.v_str(w.rc()));
        //FIXME put reachable nodes function here instead of None
        self.local_grow(&mut p, group, false, None);
        debug!("Found path {}", p.print(self.g));
        if p.in_path(v.node_id) {
            //should be covered by find_jump_ahead by this point (if possible to extend)
            return None;
        }

        //FIXME add size estimate
        p.append_general(GeneralizedLink::AMBIG(GapInfo {
            start: p.end(),
            end: v.rc(),
            gap_size: MIN_GAP_SIZE as i64}));

        Some(p.reverse_complement())
    }

    //FIXME maybe stop grow process immediately when this fails
    fn check_available(&self, node_id: usize, target_group: TrioGroup) -> bool {
        if self.in_sccs.contains(&node_id) {
            //is part of some non-trivial SCC
            return false;
        }

        if let Some(init_group) = self.assignments.group(node_id) {
            if TrioGroup::incompatible(init_group, target_group) {
                //if target group is incompatible with initial assignment (incl. ISSUE)
                return false;
            }
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

    fn local_grow(&self, path: &mut Path, group: TrioGroup, check_avail: bool,
                    //FIXME how do I typedef things like this?
                    available_vertex_f: Option<&dyn Fn(&Vertex)->bool>) -> usize {
        let mut tot_grow = 0;
        loop {
            let mut grow = self.grow_forward(path, group, true, None);
            if grow > 0 {
                debug!("Was able to extend in unambiguously assignment-aware way by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.jump_forward(path, self.find_bubble_jump_ahead(path.end(), group),
                                        group, available_vertex_f);
            if grow > 0 {
                debug!("Was able to jump across the bubble by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.jump_forward(path, self.find_small_tangle_jump_ahead(path.end(), group),
                                        group, available_vertex_f);
            if grow > 0 {
                debug!("Was able to jump across the small tangle by {grow} nodes");
                tot_grow += grow;
                continue;
            }

            assert!(grow == 0);
            break;
        }
        tot_grow
    }

    //FIXME review check_avail logic
    fn grow_forward(&self, path: &mut Path, group: TrioGroup, check_avail: bool,
                    //FIXME how do I typedef things like this?
                    available_vertex_f: Option<&dyn Fn(&Vertex)->bool>) -> usize {
        let mut v = path.end();
        let mut steps = 0;
        while let Some(l) = self.extension_helper.group_extension(v, group, available_vertex_f) {
            let w = l.end;

            if path.in_path(w.node_id)
                //FIXME relax this condition
                || self.in_sccs.contains(&w.node_id)
                || (check_avail && !self.check_available(w.node_id, group)) {
                break;
            } else {
                path.append(l);
                v = w;
                steps += 1;
            }
        }
        steps
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
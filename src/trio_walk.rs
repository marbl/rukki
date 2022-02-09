use crate::graph::*;
use crate::trio::*;
use crate::pseudo_hap::*;
use crate::graph_algos::*;
//FIXME move to common
use scc::only_or_none;
use log::debug;
use std::collections::{HashSet,HashMap};

//TODO add template parameter
pub struct HomozygousAssigner<'a> {
    g: &'a Graph,
    assignments: AssignmentStorage<'a>,
    node_len_thr: usize,
}

impl <'a> HomozygousAssigner<'a> {

    fn marking_round(&mut self) -> usize {
        //FIXME call only on the outer bubble chains
        let mut marked = 0;
        //TODO think how it should work with generalized super-bubbles
        //(probably should give a chance to extend even the node is already marked)
        for v in self.g.all_vertices() {
            //the node is long, lacks any assignment (including ISSUE) and neighborhood checks out
            if self.g.node(v.node_id).length >= self.node_len_thr
                && self.assignments.get(v.node_id).is_none()
                && self.check_homozygous_neighborhood(v) {
                marked += self.mark_vertex_and_chains(v);
            }
        }
        marked
    }

    fn mark_vertex_and_chains(&mut self, v: Vertex) -> usize {
        //hit node with existing assignment
        if !self.mark_vertex(v) {
            //already marked
            return 0;
        }
        let mut marked = 1;
        marked += self.mark_chain_ahead(v);
        marked += self.mark_chain_ahead(v.rc());
        marked
    }

    fn mark_vertex(&mut self, v: Vertex) -> bool {
        if !self.assignments.get(v.node_id).is_none() {
            //hit node with existing assignment
            return false;
        }
        self.assignments.assign(v.node_id, Assignment::<TrioGroup>{
            group: TrioGroup::HOMOZYGOUS,
            confidence: Confidence::MODERATE,
            info: String::from("HomozygousAssigner"),
        });
        true
    }

    fn mark_chain_ahead(&mut self, v: Vertex) -> usize {
        //FIXME proper parameterization
        let params = superbubble::SbSearchParams {
            max_length: 200_000,
            max_diff: 200_000,
            max_count: 1000,
        };

        let mut marked = 0;
        for bubble in superbubble::find_chain_ahead(self.g, v, &params) {
            //set to match previous logic, maybe rethink
            if !self.mark_vertex(bubble.end_vertex()) {
                break;
            }
            marked += 1;
        }
        marked
    }

    //TODO checking only one is probably enough, since iterating over all vertices
    fn check_homozygous_neighborhood(&self, v: Vertex) -> bool {
        self.check_homozygous_fork_ahead(v)
            || self.check_homozygous_fork_ahead(v.rc())
    }

    fn check_homozygous_fork_ahead(&self, v: Vertex) -> bool {
        let long_ahead = dfs::sinks_ahead(self.g, v, self.node_len_thr);
        let mut blended_group = None;
        for v_ahead in long_ahead {
            match self.assignments.group(v_ahead.node_id) {
                None | Some(TrioGroup::ISSUE) => return false,
                og => blended_group = TrioGroup::optional_blend(blended_group, og),
            };
            //looking back we should only get to v.rc()
            //TODO improve with flow ideas
            if !dfs::sinks_ahead(self.g, v_ahead.rc(), self.node_len_thr)
                .iter()
                .all(|&w| w == v.rc()) {
                return false;
            }
        }
        blended_group == Some(TrioGroup::HOMOZYGOUS)
    }
}

pub fn assign_homozygous<'a>(g: &'a Graph,
    assignments: AssignmentStorage<'a>,
    node_len_thr: usize) -> AssignmentStorage<'a> {
    let mut total_assigned = 0;
    let mut assigner = HomozygousAssigner {
        g, assignments, node_len_thr,
    };
    loop {
        debug!("Marking round");
        let marked = assigner.marking_round();
        debug!("Marked {}", marked);
        if marked == 0 {
            break;
        }
        total_assigned += marked;
    }
    debug!("Total marked {}", total_assigned);
    assigner.assignments
}

const MIN_GAP_SIZE: usize = 1000;

//TODO add template parameter
pub struct HaploSearcher<'a> {
    g: &'a Graph,
    assignments: &'a AssignmentStorage<'a>,
    long_node_threshold: usize,
    //TODO consider using same structure as for initial assignments
    used: AssignmentStorage<'a>,
    in_sccs: HashSet<usize>,
    small_tangle_index: HashMap<Vertex, scc::LocalizedTangle>,
}

impl <'a> HaploSearcher<'a> {

    pub fn new(g: &'a Graph, assignments: &'a AssignmentStorage<'a>, long_node_threshold: usize) -> HaploSearcher<'a> {
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
            used: AssignmentStorage::new(g),
            in_sccs: scc::nodes_in_sccs(g, &sccs),
            small_tangle_index,
        }
    }

    pub fn used(&self) -> &AssignmentStorage<'a> {
        &self.used
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
                let path = self.haplo_path(node_id, group);
                self.used.update_all(path.vertices().iter().map(|v| v.node_id), group);
                self.used.get_mut(path.start().node_id).unwrap().info = String::from("path_boundary");
                self.used.get_mut(path.end().node_id).unwrap().info = String::from("path_boundary");
                answer.push((path, node_id, group));
            }
        }
        answer
    }

    fn haplo_path(&self, node_id: usize, group: TrioGroup) -> Path {
        assert!(!self.incompatible_assignment(node_id, group));
        let mut path = Path::new(Vertex::forward(node_id));
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
            let mut grow = self.grow_forward(path, group, true);
            if grow > 0 {
                debug!("Was able to extend in unambiguously assignment-aware way by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.jump_forward(path, self.find_bubble_jump_ahead(path.end(), group), group);
            if grow > 0 {
                debug!("Was able to jump across the bubble by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.jump_forward(path, self.find_small_tangle_jump_ahead(path.end(), group), group);
            if grow > 0 {
                debug!("Was able to jump across the small tangle by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.jump_forward(path, self.find_jump_ahead(path.end(), group), group);
            if grow > 0 {
                debug!("Was able to jump (and stitch) ahead by {grow} nodes");
                tot_grow += grow;
                continue;
            }
            grow += self.jump_forward(path, self.find_gapped_jump_ahead(path, group), group);
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

    fn find_unbroken_alt_candidate(&self, v: Vertex) -> Option<Vertex> {
        //not necessary, but improves 'symmetry'
        assert!(self.g.vertex_length(v) >= self.long_node_threshold);

        //dead-end case
        if self.g.outgoing_edge_cnt(v) == 0 {
            let component = dfs::ShortNodeComponent::back_from_long(self.g,
                                            v, self.long_node_threshold);

            //think of maybe relaxing
            if !component.simple_boundary() {
                return None;
            }

            only_or_none(component.sinks.iter().copied().filter(|&s| s != v))
        } else if self.g.outgoing_edge_cnt(v) == 1 {
            //haplotype merge-in case
            Some(self.g.outgoing_edges(v)[0].end)
        } else {
            None
        }
    }

    //FIXME add debug prints
    //TODO very asymmetric condition :(
    fn generalized_gap(&self, v: Vertex, group: TrioGroup) -> Option<GapInfo> {
        //FIXME might be much easier to augment graph with extra 'gap' links after all!
        if self.g.vertex_length(v) < self.long_node_threshold {
            return None;
        }
        let alt = self.find_unbroken_alt_candidate(v)?;
        if self.assignments.is_definite(alt.node_id)
            && TrioGroup::incompatible(group, self.assignments.group(alt.node_id).unwrap()) {
            //debug!("Searching for component ahead from {}", self.g.v_str(alt));
            let component = dfs::ShortNodeComponent::ahead_from_long(self.g,
                                            alt, self.long_node_threshold);

            //think of maybe relaxing
            if !component.simple_boundary() {
                return None;
            }

            if component.sources.iter().all(|x| self.assignments.is_definite(x.node_id)) {
                if let Some(&w) = only_or_none(component.sources.iter()
                                .filter(|s| self.assignments.group(s.node_id).unwrap() == group)
                                .filter(|&s| self.g.incoming_edge_cnt(*s) == 0)) {
                    //dead-end case
                    //debug!("Dead-end merge-in case success");
                    return Some(GapInfo {
                        start: v,
                        end: w,
                        gap_size: std::cmp::max(self.g.vertex_length(alt) as i64
                                - self.g.vertex_length(v) as i64
                                - self.g.vertex_length(w) as i64, MIN_GAP_SIZE as i64),
                    });
                } else if component.sources.len() == 1 {
                    //haplotype merge-in case
                    assert!(component.sources.iter().next() == Some(&alt));
                    if component.sinks.iter().all(|x| self.assignments.is_definite(x.node_id)) {
                        if let Some(&w) = only_or_none(component.sinks.iter()
                                        .filter(|s| self.assignments.group(s.node_id).unwrap() == group)) {
                            //FIXME code duplication!
                            //debug!("Haplotype merge-in case success");
                            return Some(GapInfo {
                                start: v,
                                end: w,
                                gap_size: std::cmp::max(self.g.vertex_length(alt) as i64
                                        - self.g.vertex_length(v) as i64
                                        - self.g.vertex_length(w) as i64, MIN_GAP_SIZE as i64),
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
        if let Some(gap_info) = self.generalized_gap(v, group) {
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

    fn jump_forward(&self, path: &mut Path, opt_jump: Option<Path>, group: TrioGroup) -> usize {
        match opt_jump {
            None => 0,
            Some(jump) => {
                assert!(jump.len() > 1);
                assert!(path.end() == jump.start());
                //FIXME improve logging!
                if path.can_merge_in(&jump)
                    && (&(jump.vertices())[0..(jump.len() - 1)]).iter().all(|v| !self.in_sccs.contains(&v.node_id))
                    && jump.vertices().iter().all(|v| self.check_available(v.node_id, group)) {
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
            && self.assignments.contains(w.node_id)
            && !self.incompatible_assignment(w.node_id, group)
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

    fn find_jump_ahead(&self, v: Vertex, group: TrioGroup) -> Option<Path> {
        debug!("Trying to jump ahead from {}", self.g.v_str(v));
        //Currently behavior is quite conservative:
        //1. all long nodes ahead should have assignment
        //2. only one should have correct assignment
        //3. this one should have unambiguous path backward to the vertex maybe stopping one link away
        let long_ahead = dfs::sinks_ahead(self.g, v, self.long_node_threshold);

        //println!("Long ahead: {}", long_ahead.iter().map(|x| self.g.v_str(*x)).collect::<Vec<String>>().join(";"));

        //if long_ahead.iter().all(|x| self.assignments.is_definite(x.node_id)) {
        if long_ahead.iter().all(|x| self.assignments.contains(x.node_id)) {
            let potential_ext: Vec<Vertex> = long_ahead.into_iter()
                //.filter(|x| self.assignments.get(x.node_id).unwrap().group == group)
                .filter(|x| x != &v)
                .filter(|x| TrioGroup::compatible(self.assignments.group(x.node_id).unwrap(), group))
                .collect();

            debug!("Compatible extension count: {} ({})", potential_ext.len(),
                potential_ext.iter().map(|x| self.g.v_str(*x)).collect::<Vec<String>>().join(";"));

            if potential_ext.len() == 1 {
                debug!("Unique potential extension {}", self.g.v_str(potential_ext[0]));
                let mut p = Path::new(potential_ext[0].rc());
                debug!("Growing path forward from {}", self.g.v_str(potential_ext[0].rc()));
                self.grow_forward(&mut p, group, false);
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
            }
        } else {
            debug!("Not all long extensions had definite assignments");
        }

        debug!("Can't find jump");

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

    fn find_gapped_jump_ahead(&self, path: &Path, group:TrioGroup) -> Option<Path> {
        let v = path.end();
        let component = dfs::ShortNodeComponent::search_from(self.g, v, self.long_node_threshold);

        //check that sources/sinks are clearly separated and that all have assignments
        if !component.simple_boundary()
            || !component.sources.iter()
                    .chain(component.sinks.iter())
                    .all(|x| self.assignments.contains(x.node_id)) {
            return None;
        }

        let candidate_sources: Vec<Vertex> = component.sources.iter()
            .filter(|x| TrioGroup::compatible(self.assignments.group(x.node_id).unwrap(), group))
            .copied()
            .collect();

        if candidate_sources.len() != 1 || !path.vertices().contains(&candidate_sources[0]) {
            return None;
        }

        let candidate_sinks: Vec<Vertex> = component.sinks.iter()
            .filter(|x| TrioGroup::compatible(self.assignments.group(x.node_id).unwrap(), group))
            .copied()
            .collect();

        if candidate_sinks.len() != 1 {
            return None;
        }

        let w = candidate_sinks[0];

        //FIXME optimize, this info should be logged within the short node component
        if !dfs::sinks_ahead(self.g, v, self.long_node_threshold).contains(&w) {
            //w can't be reached from v
            return None;
        }

        debug!("Unique potential extension {}", self.g.v_str(w));
        let mut p = Path::new(w.rc());
        debug!("Growing path forward from {}", self.g.v_str(w.rc()));
        self.grow_forward(&mut p, group, false);
        debug!("Found path {}", p.print(self.g));
        if p.in_path(v.node_id) {
            //should be covered by find_jump_ahead by this point (if possible to extend)
            return None;
        }

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
        true
    }

    //FIXME review check_avail logic
    fn grow_forward(&self, path: &mut Path, group: TrioGroup, check_avail: bool) -> usize {
        let mut v = path.end();
        let mut steps = 0;
        while let Some(l) = self.group_extension(v, group) {
            let w = l.end;

            if path.in_path(w.node_id)
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

    fn incompatible_assignment(&self, node_id: usize, target_group: TrioGroup) -> bool {
        match self.assignments.group(node_id) {
            Some(group) => TrioGroup::incompatible(group, target_group),
            None => false,
        }
    }

    fn check_assignment(&self, node_id: usize, target_group: TrioGroup) -> bool {
        if let Some(assign) = self.assignments.get(node_id) {
            if assign.group == target_group {
                return true;
            }
        }
        false
    }

    //maybe move to graph or some GraphAlgoHelper?
    fn unambiguous_extension(&self, v: Vertex) -> Option<Link> {
        //TODO simplify?
        match self.g.outgoing_edge_cnt(v) {
            1 => Some(self.g.outgoing_edges(v)[0]),
            _ => None,
        }
    }

    //maybe move to graph or some GraphAlgoHelper?
    fn group_extension(&self, v: Vertex, group: TrioGroup) -> Option<Link> {
        if let Some(l) = self.unambiguous_extension(v) {
            if !self.incompatible_assignment(l.end.node_id, group) {
                return Some(l);
            }
        }

        let mut suitable_extension = None;
        for l in self.g.outgoing_edges(v) {
            let w = l.end;
            //currently require all extensions to be definite
            if self.assignments.is_definite(w.node_id) {
                //and exactly one having correct group assignment
                if self.assignments.group(w.node_id).unwrap() == group {
                    if suitable_extension.is_none() {
                        suitable_extension = Some(l);
                    } else {
                        return None;
                    }
                }
            } else {
                return None;
            }
        }
        suitable_extension
    }

}
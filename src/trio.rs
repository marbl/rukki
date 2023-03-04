use crate::graph::*;
use crate::graph_algos::dfs;
use crate::graph_algos::superbubble;
use log::debug;
use log::info;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Result as IOResult;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

//TODO add UNASSIGNED to display useful info for all nodes
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Eq, Ord)]
pub enum TrioGroup {
    MATERNAL,
    PATERNAL,
    HOMOZYGOUS,
    ISSUE,
}

impl TrioGroup {
    pub fn incompatible(g1: TrioGroup, g2: TrioGroup) -> bool {
        g1 == TrioGroup::ISSUE
            || g2 == TrioGroup::ISSUE
            || (g1 == TrioGroup::MATERNAL && g2 == TrioGroup::PATERNAL)
            || (g1 == TrioGroup::PATERNAL && g2 == TrioGroup::MATERNAL)
    }

    pub fn compatible(g1: TrioGroup, g2: TrioGroup) -> bool {
        !Self::incompatible(g1, g2)
    }

    pub fn is_definite(&self) -> bool {
        matches!(*self, TrioGroup::MATERNAL | TrioGroup::PATERNAL)
    }

    pub fn blend(g1: TrioGroup, g2: TrioGroup) -> TrioGroup {
        assert!(g1 != TrioGroup::ISSUE && g2 != TrioGroup::ISSUE);
        if g1 == g2 {
            g1
        } else {
            TrioGroup::HOMOZYGOUS
        }
    }

    pub fn optional_blend(og1: Option<TrioGroup>, og2: Option<TrioGroup>) -> Option<TrioGroup> {
        match og1 {
            None => og2,
            Some(g1) => match og2 {
                None => og1,
                Some(g2) => Some(Self::blend(g1, g2)),
            },
        }
    }
}

#[derive(Clone, Debug)]
pub struct Assignment {
    pub group: TrioGroup,
    pub info: String,
}

#[derive(Clone, Debug)]
pub struct TrioInfo {
    pub node_name: String,
    pub mat: usize,
    pub pat: usize,
}

impl TrioInfo {
    fn _total(&self) -> usize {
        self.mat + self.pat
    }

    fn counts_str(&self) -> String {
        format!("m{}:p{}", self.mat, self.pat)
    }
}

pub fn read_trio(path: &PathBuf) -> IOResult<Vec<TrioInfo>> {
    let mut infos = Vec::new();
    let file = File::open(path)?;
    for line in BufReader::new(file).lines() {
        let l = line?;
        let split: Vec<&str> = l.trim().split('\t').collect();
        if &split[0].to_lowercase() != "node" && &split[0].to_lowercase() != "contig" {
            let node_name = String::from(split[0]);
            let mat: usize = split[1].parse().expect("Invalid maternal count");
            let pat: usize = split[2].parse().expect("Invalid paternal count");
            infos.push(TrioInfo {
                node_name,
                mat,
                pat,
            })
        }
    }
    Ok(infos)
}

//TODO add template parameter
#[derive(Clone)]
pub struct AssignmentStorage {
    storage: HashMap<usize, Assignment>,
}

impl Default for AssignmentStorage {
    fn default() -> Self {
        Self::new()
    }
}

//TODO remove by_name methods
impl AssignmentStorage {
    pub fn new() -> AssignmentStorage {
        AssignmentStorage {
            storage: HashMap::new(),
        }
    }

    pub fn assigned(&self) -> impl Iterator<Item = usize> + '_ {
        self.storage.keys().copied()
    }

    pub fn is_definite(&self, node_id: usize) -> bool {
        if let Some(assign) = self.storage.get(&node_id) {
            if TrioGroup::is_definite(&assign.group) {
                return true;
            }
        }
        false
    }

    pub fn assign<S: Into<String>>(
        &mut self,
        node_id: usize,
        group: TrioGroup,
        info: S,
    ) -> Option<Assignment> {
        self.storage.insert(
            node_id,
            Assignment {
                group,
                info: info.into(),
            },
        )
    }

    pub fn update_group(&mut self, node_id: usize, group: TrioGroup) {
        match self.group(node_id) {
            //FIXME how to simultaneously check key and get mutable reference to stored value?
            Some(exist_group) => {
                self.storage.get_mut(&node_id).unwrap().group = TrioGroup::blend(exist_group, group)
            }
            None => {
                self.assign(node_id, group, "");
            }
        };
    }

    pub fn update_all(&mut self, iter: impl Iterator<Item = usize>, group: TrioGroup) {
        for node_id in iter {
            self.update_group(node_id, group);
        }
    }

    pub fn get(&self, node_id: usize) -> Option<&Assignment> {
        self.storage.get(&node_id)
    }

    pub fn get_mut(&mut self, node_id: usize) -> Option<&mut Assignment> {
        self.storage.get_mut(&node_id)
    }

    pub fn contains(&self, node_id: usize) -> bool {
        self.storage.contains_key(&node_id)
    }

    pub fn group(&self, node_id: usize) -> Option<TrioGroup> {
        self.storage.get(&node_id).map(|assign| assign.group)
    }
}

pub struct GroupAssignmentSettings {
    /// Minimal number of parent-specific markers required for assigning parental group to a node
    pub assign_cnt: usize,
    /// Require at least (node_length / <value>) markers within the node for parental group assignment
    pub assign_sparsity: usize,
    /// Sets minimal marker excess for assigning a parental group to <value>:1
    pub assign_ratio: f64,
    /// Sets minimal marker excess for assigning a parental group of solid nodes to <value>:1
    pub solid_ratio: f64,
    /// Minimal node length for assigning ISSUE label
    pub issue_len: usize,
    /// Minimal number of markers for assigning ISSUE label, will typically be set to a value >= assign_cnt
    pub issue_cnt: usize,
    /// Require at least (node_length / <value>) markers for assigning ISSUE label, typically set to a value >= assign_sparsity
    pub issue_sparsity: usize,
    /// Require primary marker excess BELOW <value>:1 for assigning ISSUE label. Must be <= marker_ratio
    pub issue_ratio: f64,
}

impl Default for GroupAssignmentSettings {
    fn default() -> Self {
        Self {
            assign_cnt: 10,
            assign_sparsity: 10_000,
            assign_ratio: 5.,
            solid_ratio: 5.,
            issue_len: 50_000,
            issue_cnt: 10,
            issue_sparsity: 10_000,
            issue_ratio: 5.,
        }
    }
}

pub fn assign_parental_groups(
    g: &Graph,
    trio_infos: &[TrioInfo],
    settings: &GroupAssignmentSettings,
    solid_len: usize,
) -> AssignmentStorage {
    let mut assignments = AssignmentStorage::new();

    info!("Running parental group assignment.");
    debug!("Parental group assignment settings: Minimal marker count -- {}; Minimal sparsity -- 1 in {}; Minimal ratio -- {} to 1",
            settings.assign_cnt, settings.assign_sparsity, settings.assign_ratio);
    debug!("ISSUE labeling settings: Minimal marker count -- {}; Minimal sparsity -- 1 in {}; Maximal ratio -- {} to 1",
            settings.issue_cnt, settings.issue_sparsity, settings.issue_ratio);
    assert!(settings.issue_ratio <= settings.assign_ratio);

    let assign_node_f = |x: usize, y: usize, node_len: usize| {
        assert!(x >= y);
        let tot = x + y;
        tot >= settings.assign_cnt
            && node_len <= tot * settings.assign_sparsity
            && ((x as f64) > settings.assign_ratio * (y as f64) - 1e-6
                || (node_len > solid_len && (x as f64) > settings.solid_ratio * (y as f64) - 1e-6))
    };

    let issue_node_f = |x: usize, y: usize, node_len: usize| {
        assert!(x >= y);
        let tot = x + y;
        node_len >= settings.issue_len
            && tot >= settings.issue_cnt
            && node_len <= tot * settings.issue_sparsity
            && (x as f64) < settings.issue_ratio * (y as f64) - 1e-6
    };

    for trio_info in trio_infos {
        let node_id = g.name2id(&trio_info.node_name);
        let node_len = g.node_length(node_id);
        debug!(
            "Looking at node {} (len={}), mat:pat={}",
            trio_info.node_name,
            node_len,
            trio_info.counts_str()
        );

        if issue_node_f(
            max(trio_info.mat, trio_info.pat),
            min(trio_info.mat, trio_info.pat),
            node_len,
        ) {
            debug!("Assigning ISSUE label");
            assignments.assign(node_id, TrioGroup::ISSUE, trio_info.counts_str());
        } else if assign_node_f(
            max(trio_info.mat, trio_info.pat),
            min(trio_info.mat, trio_info.pat),
            node_len,
        ) {
            if trio_info.mat >= trio_info.pat {
                debug!("Looks MATERNAL");
                assignments.assign(node_id, TrioGroup::MATERNAL, trio_info.counts_str());
            } else {
                debug!("Looks PATERNAL");
                assignments.assign(node_id, TrioGroup::PATERNAL, trio_info.counts_str());
            }
        } else {
            debug!("Failed to assign label based on marker counts");
        }
    }
    assignments
}

fn parse_group(group_str: &str) -> TrioGroup {
    match group_str {
        "MATERNAL" => TrioGroup::MATERNAL,
        "PATERNAL" => TrioGroup::PATERNAL,
        "HOMOZYGOUS" => TrioGroup::HOMOZYGOUS,
        "ISSUE" => TrioGroup::ISSUE,
        _ => panic!("Invalid group string {group_str}"),
    }
}

pub fn parse_node_assignments(
    g: &Graph,
    assignments_fn: &str,
) -> std::io::Result<AssignmentStorage> {
    let mut assignments = AssignmentStorage::new();
    for line in std::fs::read_to_string(assignments_fn)?.lines() {
        let split: Vec<&str> = line.trim().split('\t').collect();
        if &split[0].to_lowercase() != "node" && &split[0].to_lowercase() != "contig" {
            let node_name = split[0];
            let group = parse_group(split[1]);
            assignments.update_group(g.name2id(node_name), group);
        }
    }
    Ok(assignments)
}

const MAX_COMPONENT_SIZE: usize = 100;

pub struct HomozygousAssigner<'a> {
    g: &'a Graph,
    assignments: AssignmentStorage,
    trusted_len: usize,
    min_suspect_cov: Option<f64>,
    solid_len: usize,
    min_solid_cov: f64,
    max_assign_len: usize,
    complex_component_size: usize,
    considered: HashSet<usize>,
}

impl<'a> HomozygousAssigner<'a> {
    pub fn new(
        g: &'a Graph,
        assignments: AssignmentStorage,
        trusted_len: usize,
        min_suspect_cov: Option<f64>,
        solid_len: usize,
        min_solid_cov: f64,
        max_assign_len: usize,
    ) -> HomozygousAssigner<'a> {
        HomozygousAssigner {
            g,
            assignments,
            trusted_len,
            min_suspect_cov,
            solid_len,
            min_solid_cov,
            max_assign_len,
            complex_component_size: MAX_COMPONENT_SIZE,
            considered: HashSet::new(),
        }
    }

    fn can_assign(&self, node_id: usize) -> bool {
        let n = self.g.node(node_id);
        if n.length > self.max_assign_len {
            return false;
        }

        match self.assignments.group(node_id) {
            None => n.length < self.solid_len || n.coverage > self.min_solid_cov - 1e-5,
            //TODO think if we should be able to also reclassify ISSUE nodes
            Some(TrioGroup::ISSUE) => false,
            //TODO can probably be removed / asserted if only single round allowed
            Some(TrioGroup::HOMOZYGOUS) => true,
            _ => {
                n.length < self.trusted_len
                    && self.min_suspect_cov.is_some()
                    && n.coverage > self.min_suspect_cov.unwrap() - 1e-5
            }
        }
    }

    fn exclude_complicated(&mut self) {
        let mut accounted_long_starts = HashSet::new();
        for v in self.g.all_vertices() {
            if self.g.vertex_length(v) < self.trusted_len || accounted_long_starts.contains(&v) {
                continue;
            }

            let short_node_component =
                dfs::ShortNodeComponent::ahead_from_long(self.g, v, self.trusted_len);
            if short_node_component.inner.len() >= self.complex_component_size {
                for w in short_node_component.inner {
                    self.considered.insert(w.node_id);
                }
            }
            for s in &short_node_component.sources {
                accounted_long_starts.insert(*s);
            }
            for t in &short_node_component.sinks {
                accounted_long_starts.insert(t.rc());
            }
        }
    }

    pub fn run(mut self) -> AssignmentStorage {
        assert!(self.considered.is_empty());

        self.exclude_complicated();
        //FIXME call only on the outer bubble chains
        let mut marked = 0;
        //TODO think how it should work with generalized super-bubbles
        //(probably should give a chance to extend even the node is already marked)
        for v in self.g.all_vertices() {
            debug!("Considering vertex {}", self.g.v_str(v));
            if !self.considered.contains(&v.node_id)
                && self.can_assign(v.node_id)
                && self.check_homozygous_neighborhood(v)
            {
                marked += self.mark_vertex_and_chains(v);
            }
        }
        info!("Marked {}", marked);
        self.assignments
    }

    fn mark_vertex_and_chains(&mut self, v: Vertex) -> usize {
        debug!("Marking vertex {}", self.g.v_str(v));
        //hit node with existing assignment
        let mut marked = self.make_homozygous(v);
        marked += self.mark_chain_ahead(v);
        marked += self.mark_chain_ahead(v.rc());
        debug!("Done marking");
        marked
    }

    fn make_homozygous(&mut self, v: Vertex) -> usize {
        self.considered.insert(v.node_id);
        if self.can_assign(v.node_id)
            && self.assignments.group(v.node_id) != Some(TrioGroup::HOMOZYGOUS)
        {
            self.assignments
                .assign(v.node_id, TrioGroup::HOMOZYGOUS, "HomozygousAssigner");
            1
        } else {
            0
        }
    }

    fn mark_chain_ahead(&mut self, v: Vertex) -> usize {
        //FIXME proper parameterization
        let params = superbubble::SbSearchParams::unrestricted();
        let mut marked = 0;
        for bubble in superbubble::find_chain_ahead(self.g, v, &params) {
            marked += self.make_homozygous(bubble.end_vertex());
        }
        marked
    }

    //TODO checking only one is probably enough, since iterating over all vertices
    fn check_homozygous_neighborhood(&self, v: Vertex) -> bool {
        self.check_homozygous_fork_ahead(v) || self.check_homozygous_fork_ahead(v.rc())
    }

    fn check_homozygous_fork_ahead(&self, v: Vertex) -> bool {
        //trick is that v no longer has to itself be long
        let (long_ahead, mut visited_vertices) =
            dfs::sinks_ahead(self.g, v, self.trusted_len, None);
        visited_vertices.extend(&long_ahead);
        let mut blended_group = None;

        //todo maybe chack long_ahead size
        for v_ahead in &long_ahead {
            match self.assignments.group(v_ahead.node_id) {
                None | Some(TrioGroup::ISSUE) => return false,
                og => blended_group = TrioGroup::optional_blend(blended_group, og),
            };
        }

        if blended_group != Some(TrioGroup::HOMOZYGOUS) {
            return false;
        }

        //check that all incoming edges go from visited vertices
        visited_vertices.iter().all(|&x| {
            x == v
                || self
                    .g
                    .incoming_edges(x)
                    .iter()
                    .all(|&l| visited_vertices.contains(&l.start))
        })
    }
}

pub struct TangleAssignmentSettings {
    pub allow_deadend: bool,
    pub check_inner: bool,
    pub allow_reassign: bool,
}

impl Default for TangleAssignmentSettings {
    fn default() -> Self {
        Self {
            allow_deadend: false,
            check_inner: false,
            allow_reassign: true,
        }
    }
}

pub fn assign_short_node_tangles(
    g: &Graph,
    mut assignments: AssignmentStorage,
    solid_len: usize,
    settings: TangleAssignmentSettings,
) -> AssignmentStorage {
    let mut considered_boundary = HashSet::<Vertex>::new();
    for v in g.all_vertices() {
        if !considered_boundary.contains(&v) && g.vertex_length(v) >= solid_len {
            let comp = dfs::ShortNodeComponent::ahead_from_long(g, v, solid_len);

            for s in comp.sources.iter() {
                considered_boundary.insert(*s);
            }
            for s in comp.sinks.iter() {
                considered_boundary.insert(s.rc());
            }

            if !settings.allow_deadend && comp.has_deadends {
                continue;
            }

            if !assignments.is_definite(v.node_id) {
                continue;
            }
            let group = assignments.group(v.node_id).unwrap();

            if comp
                .sources
                .iter()
                .chain(comp.sinks.iter())
                .any(|&x| Some(group) != assignments.group(x.node_id))
            {
                continue;
            }

            if settings.check_inner
                && comp
                    .inner
                    .iter()
                    .any(|&x| group != assignments.group(x.node_id).unwrap_or(group))
            {
                continue;
            }

            for w in comp.inner.iter() {
                match assignments.group(w.node_id) {
                    None => {
                        assignments.assign(w.node_id, group, "TangleAssignment");
                    }
                    Some(g) if g == group => {}
                    _ => {
                        if settings.allow_reassign {
                            assignments.assign(w.node_id, group, "TangleReAssignment");
                        }
                    }
                }
            }
        }
    }
    assignments
}

#[cfg(test)]
mod tests {
    use crate::graph::*;
    use crate::trio;
    use std::fs;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn homozygous_fork_test() {
        init();

        let graph_fn = "tests/test_graphs/test1.gfa";
        let assignments_fn = "tests/test_graphs/test1.ann.csv";
        let g = Graph::read(&fs::read_to_string(graph_fn).unwrap());
        let assignments = trio::parse_node_assignments(&g, assignments_fn).unwrap();

        let assigner =
            trio::HomozygousAssigner::new(&g, assignments, 100_000, None, 500_000, 1.5, usize::MAX);
        assert!(assigner.check_homozygous_fork_ahead(Vertex::forward(g.name2id("utig4-1237"))));
        assert!(assigner.check_homozygous_fork_ahead(Vertex::reverse(g.name2id("utig4-1237"))));
        assert!(!assigner.check_homozygous_fork_ahead(Vertex::forward(g.name2id("utig4-1554"))));
        assert!(!assigner.check_homozygous_fork_ahead(Vertex::reverse(g.name2id("utig4-1554"))));
    }
}

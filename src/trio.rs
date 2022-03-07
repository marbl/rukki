use log::debug;
use std::collections::{HashMap,HashSet};
use log::info;
use crate::graph::*;
use crate::graph_algos::dfs;
use crate::graph_algos::superbubble;
use std::cmp::{min, max};

//#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
//pub enum Confidence {
//    INCONCLUSIVE,
//    LOW,
//    MODERATE,
//    HIGH,
//}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Eq, Ord)]
pub enum TrioGroup {
    MATERNAL,
    PATERNAL,
    HOMOZYGOUS,
    ISSUE,
}

impl TrioGroup {
    pub fn incompatible(g1: TrioGroup, g2: TrioGroup) -> bool {
        g1 == TrioGroup::ISSUE || g2 == TrioGroup::ISSUE ||
            (g1 == TrioGroup::MATERNAL && g2 == TrioGroup::PATERNAL) ||
            (g1 == TrioGroup::PATERNAL && g2 == TrioGroup::MATERNAL)
    }

    pub fn compatible(g1: TrioGroup, g2: TrioGroup) -> bool {
        !Self::incompatible(g1, g2)
    }

    pub fn is_definite(&self) -> bool {
        match *self {
            TrioGroup::MATERNAL | TrioGroup::PATERNAL => true,
            _ => false,
        }
    }

    pub fn blend(g1: TrioGroup, g2: TrioGroup) -> TrioGroup {
        assert!(g1 != TrioGroup::ISSUE && g2 != TrioGroup::ISSUE);
        if g1 == g2 {
            g1
        } else {
            return TrioGroup::HOMOZYGOUS;
        }
    }

    pub fn optional_blend(og1: Option<TrioGroup>, og2: Option<TrioGroup>) -> Option<TrioGroup> {
        match og1 {
            None => og2,
            Some(g1) => match og2 {
                None => og1,
                Some(g2) => Some(Self::blend(g1, g2)),
            }
        }
    }
}

#[derive(Clone, Debug)]
pub struct Assignment<T>
{
    pub group: T,
    //pub confidence: Confidence,
    pub info: String,
}

//impl <T> Assignment<T> {
//    fn new_basic(group: T) -> Assignment<T> {
//        Assignment {
//            group,
//            //confidence: Confidence::MODERATE,
//            info: String::new(),
//        }
//    }
//}

#[derive(Clone, Debug)]
pub struct TrioInfo {
    node_name: String,
    mat: usize,
    pat: usize,
}

impl TrioInfo {
    fn _total(&self) -> usize {
        self.mat + self.pat
    }

    fn counts_str(&self) -> String {
        format!("m{}:p{}", self.mat, self.pat)
    }
}

pub fn read_trio(trio_str: &str) -> Vec<TrioInfo> {
    let mut infos = Vec::new();
    for line in trio_str.lines() {
        let split: Vec<&str> = line.trim().split('\t').collect();
        if &split[0].to_lowercase() != "node" && &split[0].to_lowercase() != "contig" {
            let node_name = String::from(split[0]);
            let mat: usize = split[1].parse().expect("Invalid maternal count");
            let pat: usize = split[2].parse().expect("Invalid paternal count");
            infos.push(TrioInfo{node_name, mat, pat})
        }
    }
    infos
}

//TODO add template parameter
#[derive(Clone)]
pub struct AssignmentStorage {
    storage: HashMap<usize, Assignment<TrioGroup>>,
}

//TODO remove by_name methods
impl <'a> AssignmentStorage {
    pub fn new() -> AssignmentStorage {
        AssignmentStorage {
            storage: HashMap::new(),
        }
    }

    pub fn assigned(&self) -> impl Iterator<Item=usize> + '_ {
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

    pub fn assign(&mut self, node_id: usize, group: TrioGroup, info: String) -> Option<Assignment::<TrioGroup>> {
        self.storage.insert(node_id, Assignment::<TrioGroup>{group, info})
    }

    pub fn update_group(&mut self, node_id: usize, group: TrioGroup) {
        match self.group(node_id) {
            //FIXME how to simultaneously check key and get mutable reference to stored value?
            Some(exist_group) => {
                self.storage.get_mut(&node_id).unwrap().group = TrioGroup::blend(exist_group, group)
            }
            None => {
                self.assign(node_id, group, String::new());
            }
        };
    }

    pub fn update_all(&mut self, iter: impl Iterator<Item=usize>, group: TrioGroup) {
        for node_id in iter {
            self.update_group(node_id, group);
        }
    }

    //pub fn assign_by_name(&mut self, node_name: &str, assignment: Assignment<TrioGroup>) {
    //    self.assign(self.g.name2id(node_name), assignment);
    //}

    pub fn get(&self, node_id: usize) -> Option<&Assignment<TrioGroup>> {
        self.storage.get(&node_id)
    }

    pub fn get_mut(&mut self, node_id: usize) -> Option<&mut Assignment<TrioGroup>> {
        self.storage.get_mut(&node_id)
    }

    pub fn contains(&self, node_id: usize) -> bool {
        self.storage.contains_key(&node_id)
    }

    pub fn group(&self, node_id: usize) -> Option<TrioGroup> {
        if let Some(assign) = self.storage.get(&node_id) {
            Some(assign.group)
        } else {
            None
        }
    }

}

pub fn assign_parental_groups(g: &Graph, trio_infos: &[TrioInfo],
        assign_cnt: usize, assign_sparsity: usize, assign_ratio: f64,
        issue_cnt: usize, issue_sparsity: usize, issue_ratio: f64) -> AssignmentStorage {
    let mut assignments = AssignmentStorage::new();

    info!("Running parental group assignment.");
    debug!("Parental group assignment settings. Minimal marker count: {assign_cnt}; Minimal sparsity: 1 in {assign_sparsity}; Minimal ratio: {assign_ratio} to 1");
    debug!("ISSUE labeling settings. Minimal marker count: {issue_cnt}; Minimal sparsity: 1 in {issue_sparsity}; Maximal ratio: {issue_ratio} to 1");
    assert!(issue_ratio <= assign_ratio);

    let assign_node_f = |x: usize, y: usize, node_len: usize| {
        assert!(x >= y);
        let tot = x + y;
        tot >= assign_cnt && node_len <= tot * assign_sparsity
            && (x as f64) > assign_ratio * (y as f64) - 1e-6
    };

    let issue_node_f = |x: usize, y: usize, node_len: usize| {
        assert!(x >= y);
        let tot = x + y;
        tot >= issue_cnt && node_len <= tot * issue_sparsity
            && (x as f64) < issue_ratio * (y as f64) - 1e-6
    };

    for trio_info in trio_infos {
        let node_id = g.name2id(&trio_info.node_name);
        let node_len = g.node_length(node_id);
        debug!("Looking at node {} (len={}), mat:pat={}",
            trio_info.node_name, node_len, trio_info.counts_str());

        if issue_node_f(max(trio_info.mat, trio_info.pat),
            min(trio_info.mat, trio_info.pat),
            node_len) {
            debug!("Assigning ISSUE label");
            assignments.assign(node_id, TrioGroup::ISSUE, trio_info.counts_str());
        } else if assign_node_f(max(trio_info.mat, trio_info.pat),
            min(trio_info.mat, trio_info.pat),
            node_len) {
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

//pub fn old_assign_parental_groups<'a>(g: &'a Graph, trio_infos: &[TrioInfo],
//        assign_cnt: usize, assign_sparsity: usize, assign_ratio: f32) -> AssignmentStorage<'a> {
//    let mut assignments = AssignmentStorage::new(g);
//    let moderate_cnt_thr = assign_cnt * 10;
//    let high_cnt_thr =  moderate_cnt_thr * 10;
//
//    debug!("Running parental group assignment with marker confidence thresholds: high={high_cnt_thr},
//medium={moderate_cnt_thr}, low={assign_cnt}; Maximal sparsity: 1:{assign_sparsity}; Ratio: {assign_ratio}");
//    assert!(high_cnt_thr as f32 / assign_cnt as f32 > assign_ratio);
//
//    let issue_confidence = |x: usize| {
//        if x > high_cnt_thr {
//            Confidence::HIGH
//        } else if x > moderate_cnt_thr {
//            Confidence::MODERATE
//        } else if x > assign_cnt {
//            Confidence::LOW
//        } else {
//            Confidence::INCONCLUSIVE
//        }
//    };
//
//    let excess_confidence = |x: usize, y: usize| {
//        assert!(x >= y);
//        if x > high_cnt_thr && y < assign_cnt {
//            return Confidence::HIGH;
//        }
//        if y == 0 || (x as f32 / y as f32) >= assign_ratio {
//            if x > moderate_cnt_thr {
//                return Confidence::MODERATE;
//            } else if x > assign_cnt {
//                return Confidence::LOW;
//            }
//        }
//        Confidence::INCONCLUSIVE
//    };
//
//    let classify_cnts = |x: usize, y: usize| {
//        assert!(x >= y);
//        match excess_confidence(x, y) {
//            Confidence::INCONCLUSIVE => (None, issue_confidence(x)),
//            conf => (Some(conf), Confidence::INCONCLUSIVE),
//        }
//    };
//
//    for trio_info in trio_infos {
//        let node_id = g.name2id(&trio_info.node_name);
//        let node_len = g.node_length(node_id);
//        debug!("Looking at node {} (len={}), mat:pat={}",
//            trio_info.node_name, node_len, trio_info.counts_str());
//
//        //TODO maybe take max?
//        if trio_info.total() < moderate_cnt_thr && node_len > trio_info.total() * assign_sparsity {
//            debug!("Marker density lower than 1 / {assign_sparsity}");
//            continue;
//        }
//
//        //FIXME refactor!
//        if trio_info.mat >= trio_info.pat {
//            debug!("Looks more maternal");
//            match classify_cnts(trio_info.mat, trio_info.pat) {
//                (None, Confidence::INCONCLUSIVE) => {},
//                (None, issue_conf) => assignments.assign(node_id,
//                    Assignment::<TrioGroup>{
//                        group: TrioGroup::ISSUE,
//                        confidence: issue_conf,
//                        //info: String::from("Marker issue")}),
//                        info: trio_info.counts_str(),
//                    }),
//                (Some(conf), _) => assignments.assign(node_id,
//                    Assignment::<TrioGroup>{
//                        group: TrioGroup::MATERNAL,
//                        confidence: conf,
//                        //info: String::from("Maternal marker excess")}),
//                        info: trio_info.counts_str(),
//                    }),
//            }
//        } else {
//            debug!("Looks more paternal");
//            match classify_cnts(trio_info.pat, trio_info.mat) {
//                (None, Confidence::INCONCLUSIVE) => {},
//                (None, issue_conf) => assignments.assign(node_id,
//                    Assignment::<TrioGroup>{
//                        group: TrioGroup::ISSUE,
//                        //confidence: issue_conf,
//                        //info: String::from("Marker issue")}),
//                        info: trio_info.counts_str(),
//                    }),
//                (Some(conf), _) => assignments.assign(node_id,
//                    Assignment::<TrioGroup>{
//                        group: TrioGroup::PATERNAL,
//                        //confidence: conf,
//                        //info: String::from("Paternal marker excess")}),
//                        info: trio_info.counts_str(),
//                    }),
//            }
//        }
//    }
//    assignments
//}

fn parse_group(group_str: &str) -> TrioGroup {
    match group_str {
        "MATERNAL" => TrioGroup::MATERNAL,
        "PATERNAL" => TrioGroup::PATERNAL,
        "HOMOZYGOUS" => TrioGroup::HOMOZYGOUS,
        "ISSUE" => TrioGroup::ISSUE,
        _ => panic!("Invalid group string {group_str}"),
    }
}

pub fn parse_read_assignments(g: &Graph, assignments_fn: &str)
-> std::io::Result<AssignmentStorage> {
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

//TODO add template parameter
pub struct HomozygousAssigner<'a> {
    g: &'a Graph,
    assignments: AssignmentStorage,
    node_len_thr: usize,
    considered: HashSet<usize>,
}

impl <'a> HomozygousAssigner<'a> {

    pub fn new(g: &'a Graph, assignments: AssignmentStorage, node_len_thr: usize)
    -> HomozygousAssigner<'a> {
        HomozygousAssigner {
            g, assignments, node_len_thr,
            considered: HashSet::new(),
        }
    }

    fn exclude_complicated(&mut self, max_component_size: usize) {
        let mut accounted_long_starts = HashSet::new();
        for v in self.g.all_vertices() {
            if self.g.vertex_length(v) < self.node_len_thr
                || accounted_long_starts.contains(&v) {
                continue;
            }

            let short_node_component = dfs::ShortNodeComponent::ahead_from_long(self.g, v, self.node_len_thr);
            if short_node_component.inner.len() > max_component_size {
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

    fn marking_round(&mut self) -> usize {
        self.considered.clear();
        self.exclude_complicated(100);
        //FIXME call only on the outer bubble chains
        let mut marked = 0;
        //TODO think how it should work with generalized super-bubbles
        //(probably should give a chance to extend even the node is already marked)
        for v in self.g.all_vertices() {
            debug!("Considering vertex {}", self.g.v_str(v));
            if !self.considered.contains(&v.node_id)
                && (self.assignments.get(v.node_id).is_none()
                    || self.assignments.group(v.node_id) == Some(TrioGroup::HOMOZYGOUS))
                && self.check_homozygous_neighborhood(v) {
                marked += self.mark_vertex_and_chains(v);
            }
        }
        marked
    }

    fn mark_vertex_and_chains(&mut self, v: Vertex) -> usize {
        debug!("Marking vertex {}", self.g.v_str(v));
        //hit node with existing assignment
        let mut marked = self.process_homozygous_vertex(v);
        marked += self.mark_chain_ahead(v);
        marked += self.mark_chain_ahead(v.rc());
        debug!("Done marking");
        marked
    }

    fn process_homozygous_vertex(&mut self, v: Vertex) -> usize {
        self.considered.insert(v.node_id);
        if self.assignments.get(v.node_id).is_none() {
            self.assignments.assign(v.node_id,
                TrioGroup::HOMOZYGOUS,
                String::from("HomozygousAssigner"));
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
            marked += self.process_homozygous_vertex(bubble.end_vertex());
        }
        marked
    }

    //TODO checking only one is probably enough, since iterating over all vertices
    fn check_homozygous_neighborhood(&self, v: Vertex) -> bool {
        self.check_homozygous_fork_ahead(v)
            || self.check_homozygous_fork_ahead(v.rc())
    }

    fn check_homozygous_fork_ahead(&self, v: Vertex) -> bool {
        //trick is that v no longer has to itself be long
        let (long_ahead, mut visited_vertices) = dfs::sinks_ahead(self.g, v, self.node_len_thr);
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

        for &v_ahead in &long_ahead {
            let mut dfs = dfs::DFS::new_reverse(self.g);
            dfs.set_max_node_len(self.node_len_thr);
            dfs.extend_blocked(std::iter::once(v));
            dfs.run_from(v_ahead);

            //looking back we should only get to v.rc()
            //TODO improve with flow ideas
            //TODO improve performance if necessary
            if dfs.blocked().iter().chain(dfs.boundary().iter())
                .any(|x| !visited_vertices.contains(x)) {
                return false;
            }
        }
        true
    }
}

pub fn assign_homozygous(g: &Graph,
    assignments: AssignmentStorage,
    node_len_thr: usize) -> AssignmentStorage {
    info!("Marking homozygous nodes");
    let mut total_assigned = 0;
    let mut assigner = HomozygousAssigner::new(g, assignments, node_len_thr);
    loop {
        info!("Marking round");
        let marked = assigner.marking_round();
        info!("Marked {}", marked);
        if marked == 0 {
            break;
        }
        total_assigned += marked;
    }
    info!("Total marked {}", total_assigned);
    assigner.assignments
}

////FIXME two thresholds?
//pub fn assign_small_bubbles<'a>(g: &'a Graph,
//    mut assignments: AssignmentStorage<'a>,
//    node_len_thr: usize) -> AssignmentStorage<'a> {
//    let end_cov = |l: &Link| {g.node(l.end.node_id).coverage};
//
//    for v in g.all_vertices() {
//        if g.vertex_length(v) < node_len_thr || !assignments.is_definite(v.node_id) {
//            continue;
//        }
//        if let Some(w) = superbubble::trivial_bubble_end(g, v) {
//            if g.outgoing_edges(v).iter().all(|&l| g.vertex_length(l.end) < node_len_thr
//                                                    && !assignments.contains(l.end.node_id))
//                && w.node_id != v.node_id {
//                let group = assignments.group(v.node_id).unwrap();
//                for l in g.outgoing_edges(v) {
//                    assignments.assign(l.end.node_id,
//                            Assignment::<TrioGroup>{
//                                        group,
//                                        confidence: Confidence::MODERATE,
//                                        info: String::from("SmallBubble"),
//                    });
//                }
//                //assignments.assign(g.outgoing_edges(v).iter()
//                //                            .max_by(|a, b| end_cov(a).partial_cmp(&end_cov(b)).unwrap())
//                //                            .unwrap().end.node_id,
//                //        Assignment::<TrioGroup>{
//                //                    group,
//                //                    confidence: Confidence::MODERATE,
//                //                    info: String::from("SmallBubble"),
//                //});
//            }
//        }
//    }
//    assignments
//}

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
        let assignments = trio::parse_read_assignments(&g, assignments_fn).unwrap();

        let assigner = trio::HomozygousAssigner::new(&g, assignments, 100_000);
        assert!(assigner.check_homozygous_fork_ahead(Vertex::forward(g.name2id("utig4-1237"))));
        assert!(assigner.check_homozygous_fork_ahead(Vertex::reverse(g.name2id("utig4-1237"))));
        assert!(!assigner.check_homozygous_fork_ahead(Vertex::forward(g.name2id("utig4-1554"))));
        assert!(!assigner.check_homozygous_fork_ahead(Vertex::reverse(g.name2id("utig4-1554"))));
    }

}
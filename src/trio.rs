use log::debug;
use std::collections::HashMap;

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub enum Confidence {
    INCONCLUSIVE,
    LOW,
    MODERATE,
    HIGH,
}

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
    pub confidence: Confidence,
    pub info: String,
}

impl <T> Assignment<T> {
    fn new_basic(group: T) -> Assignment<T> {
        Assignment {
            group,
            confidence: Confidence::MODERATE,
            info: String::new(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct TrioInfo {
    node_name: String,
    mat: usize,
    pat: usize,
}

impl TrioInfo {
    fn total(&self) -> usize {
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

use crate::graph::Graph;

//TODO add template parameter
#[derive(Clone)]
pub struct AssignmentStorage<'a> {
    storage: HashMap<usize, Assignment<TrioGroup>>,
    g: &'a Graph,
}

//TODO remove by_name methods
impl <'a> AssignmentStorage<'a> {
    pub fn new(g: &'a Graph) -> AssignmentStorage<'a> {
        AssignmentStorage {
            storage: HashMap::new(),
            g,
        }
    }

    pub fn is_definite(&self, node_id: usize) -> bool {
        if let Some(assign) = self.storage.get(&node_id) {
            if TrioGroup::is_definite(&assign.group) {
                return true;
            }
        }
        false
    }

    pub fn assign(&mut self, node_id: usize, assignment: Assignment<TrioGroup>) {
        self.storage.insert(node_id, assignment);
    }

    pub fn update_group(&mut self, node_id: usize, group: TrioGroup) {
        match self.group(node_id) {
            //FIXME how to simultaneously check key and get mutable reference to stored value?
            Some(exist_group) => {
                self.storage.get_mut(&node_id).unwrap().group = TrioGroup::blend(exist_group, group)
            }
            None => {
                self.assign(node_id, Assignment::<TrioGroup>::new_basic(group));
            }
        };
    }

    pub fn update_all(&mut self, iter: impl Iterator<Item=usize>, group: TrioGroup) {
        for node_id in iter {
            self.update_group(node_id, group);
        }
    }

    pub fn assign_by_name(&mut self, node_name: &str, assignment: Assignment<TrioGroup>) {
        self.assign(self.g.name2id(node_name), assignment);
    }

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

pub fn assign_parental_groups<'a>(g: &'a Graph, trio_infos: &[TrioInfo], low_cnt_thr : usize, ratio_thr : f32) -> AssignmentStorage<'a> {
    let mut assignments = AssignmentStorage::new(g);
    let min_marker_inv_density = 10_000;
    let moderate_cnt_thr = low_cnt_thr * 10;
    let high_cnt_thr =  moderate_cnt_thr * 10;

    debug!("Running with marker parameters confidence: high={} medium={}, low={}; Ratio: {}",
            high_cnt_thr, moderate_cnt_thr, low_cnt_thr, ratio_thr);
    assert!(high_cnt_thr as f32 / low_cnt_thr as f32 > ratio_thr);

    let issue_confidence = |x: usize| {
        if x > high_cnt_thr {
            Confidence::HIGH
        } else if x > moderate_cnt_thr {
            Confidence::MODERATE
        } else if x > low_cnt_thr {
            Confidence::LOW
        } else {
            Confidence::INCONCLUSIVE
        }
    };

    let excess_confidence = |x: usize, y: usize| {
        assert!(x >= y);
        if x > high_cnt_thr && y < low_cnt_thr {
            return Confidence::HIGH;
        }
        if y == 0 || (x as f32 / y as f32) >= ratio_thr {
            if x > moderate_cnt_thr {
                return Confidence::MODERATE;
            } else if x > low_cnt_thr {
                return Confidence::LOW;
            }
        }
        Confidence::INCONCLUSIVE
    };

    let classify_cnts = |x: usize, y: usize| {
        assert!(x >= y);
        match excess_confidence(x, y) {
            Confidence::INCONCLUSIVE => (None, issue_confidence(x)),
            conf => (Some(conf), Confidence::INCONCLUSIVE),
        }
    };

    for trio_info in trio_infos {
        let node_len = g.node_by_name(&trio_info.node_name).length;
        debug!("Looking at node {} (len={}), mat:pat={}",
            trio_info.node_name, node_len, trio_info.counts_str());

        //TODO maybe take max?
        if trio_info.total() < moderate_cnt_thr && node_len > trio_info.total() * min_marker_inv_density {
            //FIXME continue?!!!
            debug!("Too few markers")
        }

        if trio_info.mat >= trio_info.pat {
            debug!("Looks more maternal");
            match classify_cnts(trio_info.mat, trio_info.pat) {
                (None, Confidence::INCONCLUSIVE) => {},
                (None, issue_conf) => assignments.assign_by_name(&trio_info.node_name,
                    Assignment::<TrioGroup>{
                        group: TrioGroup::ISSUE,
                        confidence: issue_conf,
                        //info: String::from("Marker issue")}),
                        info: trio_info.counts_str(),
                    }),
                (Some(conf), _) => assignments.assign_by_name(&trio_info.node_name,
                    Assignment::<TrioGroup>{
                        group: TrioGroup::MATERNAL,
                        confidence: conf,
                        //info: String::from("Maternal marker excess")}),
                        info: trio_info.counts_str(),
                    }),
            }
        } else {
            debug!("Looks more paternal");
            match classify_cnts(trio_info.pat, trio_info.mat) {
                (None, Confidence::INCONCLUSIVE) => {},
                (None, issue_conf) => assignments.assign_by_name(&trio_info.node_name,
                    Assignment::<TrioGroup>{
                        group: TrioGroup::ISSUE,
                        confidence: issue_conf,
                        //info: String::from("Marker issue")}),
                        info: trio_info.counts_str(),
                    }),
                (Some(conf), _) => assignments.assign_by_name(&trio_info.node_name,
                    Assignment::<TrioGroup>{
                        group: TrioGroup::PATERNAL,
                        confidence: conf,
                        //info: String::from("Paternal marker excess")}),
                        info: trio_info.counts_str(),
                    }),
            }
        }
    }
    assignments
}

fn parse_group(group_str: &str) -> TrioGroup {
    match group_str {
        "MATERNAL" => TrioGroup::MATERNAL,
        "PATERNAL" => TrioGroup::PATERNAL,
        "HOMOZYGOUS" => TrioGroup::HOMOZYGOUS,
        _ => panic!("Invalid group string {group_str}"),
    }
}

pub fn parse_read_assignments<'a>(g: &'a Graph, assignments_fn: &str, load_homozygous: bool)
-> std::io::Result<AssignmentStorage<'a>> {
    let mut assignments = AssignmentStorage::new(g);
    for line in std::fs::read_to_string(assignments_fn)?.lines() {
        let split: Vec<&str> = line.trim().split('\t').collect();
        if &split[0].to_lowercase() != "node" && &split[0].to_lowercase() != "contig" {
            let node_name = split[0];
            let group = parse_group(split[1]);
            if load_homozygous || group != TrioGroup::HOMOZYGOUS {
                assignments.update_group(g.name2id(node_name), group);
            }
        }
    }
    Ok(assignments)
}

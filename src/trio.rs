use std::collections::HashMap;

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub enum Confidence {
    INCONCLUSIVE,
    LOW,
    MODERATE,
    HIGH,
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
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
}

#[derive(Clone, Debug)]
pub struct Assignment<T>
where T: Clone
{
    pub group: T,
    pub confidence: Confidence,
    pub info: String,
}

#[derive(Clone, Debug)]
pub struct TrioInfo {
    node_name: String,
    mat: u64,
    pat: u64,
}

impl TrioInfo {
    fn counts_str(&self) -> String {
        format!("{}:{}", self.mat, self.pat)
    }
}

pub fn read_trio(trio_str: &str) -> Vec<TrioInfo> {
    let mut infos = Vec::new();
    for line in trio_str.lines() {
        let split: Vec<&str> = line.trim().split('\t').collect();
        if &split[0].to_lowercase() != "node" && &split[0].to_lowercase() != "contig" {
            let node_name = String::from(split[0]);
            let mat: u64 = split[1].parse().expect("Invalid maternal count");
            let pat: u64 = split[2].parse().expect("Invalid paternal count");
            infos.push(TrioInfo{node_name, mat, pat})
        }
    }
    infos
}

use crate::graph::Graph;

//TODO add template parameter
pub struct AssignmentStorage<'a> {
    storage: HashMap<usize, Assignment<TrioGroup>>,
    g: &'a Graph,
}

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

    fn assign(&mut self, node_id: usize, assignment: Assignment<TrioGroup>) {
        self.storage.insert(node_id, assignment);
    }

    fn assign_by_name(&mut self, node_name: &str, assignment: Assignment<TrioGroup>) {
        self.assign(self.g.name2id(node_name), assignment);
    }

    pub fn get(&self, node_id: usize) -> Option<&Assignment<TrioGroup>> {
        self.storage.get(&node_id)
    }

    pub fn get_by_name(&self, node_name: &str) -> Option<&Assignment<TrioGroup>> {
        self.get(self.g.name2id(node_name))
    }
}

pub fn assign_parental_groups<'a>(g: &'a Graph, trio_infos: &[TrioInfo]) -> AssignmentStorage<'a> {
    let mut assignments = AssignmentStorage::new(g);
    let high_cnt_thr = 1000;
    let moderate_cnt_thr = 300;
    let low_cnt_thr = 100;
    let ratio_thr = 5.0;
    assert!(high_cnt_thr as f32 / low_cnt_thr as f32 > ratio_thr);

    let issue_confidence = |x: u64| {
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

    let excess_confidence = |x: u64, y: u64| {
        assert!(x >= y);
        if x > high_cnt_thr && y < low_cnt_thr {
            return Confidence::HIGH;
        }
        if (x as f32 / y as f32) >= ratio_thr {
            if x > moderate_cnt_thr {
                return Confidence::MODERATE;
            } else if x > low_cnt_thr {
                return Confidence::LOW;
            }
        }
        Confidence::INCONCLUSIVE
    };

    let classify_cnts = |x: u64, y: u64| {
        assert!(x >= y);
        match excess_confidence(x, y) {
            Confidence::INCONCLUSIVE => (None, issue_confidence(x)),
            conf => (Some(conf), Confidence::INCONCLUSIVE),
        }
    };

    for trio_info in trio_infos {
        if trio_info.mat >= trio_info.pat {
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
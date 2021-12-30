use std::collections::HashMap;
use std::collections::HashSet;

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
    fn incompatible(g1: TrioGroup, g2: TrioGroup) -> bool {
        g1 == TrioGroup::ISSUE || g2 == TrioGroup::ISSUE ||
            (g1 == TrioGroup::MATERNAL && g2 == TrioGroup::PATERNAL) ||
            (g1 == TrioGroup::PATERNAL && g2 == TrioGroup::MATERNAL)
    }

    fn is_definite(&self) -> bool {
        match *self {
            TrioGroup::MATERNAL | TrioGroup::PATERNAL => true,
            _ => false,
        }
    }

    fn blend(g1: TrioGroup, g2: TrioGroup) -> TrioGroup {
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

    fn is_definite(&self, node_id: usize) -> bool {
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
    let moderate_cnt_thr = 500;
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

use crate::graph::*;

//TODO add template parameter
pub struct HaploPathSearcher<'a> {
    g: &'a Graph,
    assignments: &'a AssignmentStorage<'a>,
    long_node_threshold: usize,
    used: HashMap<usize, TrioGroup>,
}

impl <'a> HaploPathSearcher<'a> {
    pub fn new(g: &'a Graph, assignments: &'a AssignmentStorage<'a>, long_node_threshold: usize) -> HaploPathSearcher<'a> {
        HaploPathSearcher{ g, assignments, long_node_threshold, used: HashMap::new() }
    }

    fn update_used(&mut self, path: &HaploPath, group: TrioGroup) {
        for v in path.vertices() {
            let blended = match self.used.get(&v.node_id) {
                Some(exist_group) => TrioGroup::blend(*exist_group, group),
                None => group,
            };
            self.used.insert(v.node_id, blended);
        }
    }

    //TODO maybe use single length threshold?
    pub fn find_all(&mut self) -> Vec<(HaploPath, TrioGroup)> {
        let mut answer = Vec::new();
        for (node_id, node) in self.g.all_nodes().enumerate() {
            if self.used.contains_key(&node_id) {
                continue;
            }
            if node.length >= self.long_node_threshold && self.assignments.is_definite(node_id) {
                let group = self.assignments.get(node_id).unwrap().group;
                let path = self.haplo_path(node_id, group);
                self.update_used(&path, group);
                answer.push((path, group));
            }
        }
        answer
    }

    fn haplo_path(&self, node_id: usize, group: TrioGroup) -> HaploPath {
        assert!(!self.incompatible_assignment(node_id, group));
        let mut path = HaploPath::new(Vertex{node_id, direction: Direction::FORWARD});
        self.grow_forward(&mut path, group);
        path = path.reverse_complement();
        self.grow_forward(&mut path, group);
        path.reverse_complement()
    }

    //TODO maybe consume when grow?
    fn grow_forward(&self, path: &mut HaploPath, group: TrioGroup) -> usize {
        let mut tot_grow = self.unambig_grow_forward(path, group);
        //println!("Trying to jump ahead from {}", self.g.v_str(*path.end().unwrap()));
        while let Some(jump) = self.jump_ahead(*path.end().unwrap(), group) {
            //println!("Successful jump to {}", self.g.v_str(*jump.end().unwrap()));
            assert!(path.end().unwrap() == jump.start().unwrap());
            if path.can_merge_in(&jump) {
                tot_grow += jump.len() - 1;
                path.merge_in(jump);
            }
            tot_grow += self.unambig_grow_forward(path, group);
            //println!("Trying to jump ahead from {}", self.g.v_str(*path.end().unwrap()));
        }
        tot_grow
    }

    fn inner_dfs(&self, v: Vertex, visited: &mut HashSet<Vertex>, long_ext: &mut Vec<Vertex>) {
        visited.insert(v);
        //if only one vertex is visited then it means we just started
        if visited.len() > 1 && self.g.node(v.node_id).length >= self.long_node_threshold {
            long_ext.push(v);
        } else {
            for l in self.g.outgoing_edges(v) {
                let w = l.end;
                if !visited.contains(&w) {
                    self.inner_dfs(w, visited, long_ext);
                }
            }
        }
    }

    fn bounded_dfs(&self, v: Vertex) -> Vec<Vertex> {
        let mut visited = HashSet::new();
        let mut long_ext = Vec::new();
        self.inner_dfs(v, &mut visited, &mut long_ext);
        long_ext
    }

    fn try_link(&self, mut path: HaploPath, v: Vertex) -> HaploPath {
        let end = path.end().unwrap();
        for l in self.g.outgoing_edges(*end) {
            if l.end == v {
                path.append(l);
                break;
            }
        }
        path
    }


    fn try_link_with_vertex(&self, mut path: HaploPath, v: Vertex, group: TrioGroup) -> HaploPath {
        let end = path.end().unwrap();
        for l in self.g.outgoing_edges(*end) {
            let w = l.end;
            if path.in_path(w.node_id)
                //TODO think if checking length is necessary here
                || self.g.node(w.node_id).length >= self.long_node_threshold
                || self.incompatible_assignment(w.node_id, group)
                || self.g.incoming_edge_cnt(w) != 1
                || self.g.outgoing_edge_cnt(w) != 1 {
                //FIXME think if we should check coverage too
                continue;
            }
            if let Some(l2) = self.g.connector(w, v) {
                path.append(l);
                path.append(l2);
                break;
            }
        }
        path
    }

    fn jump_ahead(&self, v: Vertex, group: TrioGroup) -> Option<HaploPath> {
        //Currently behavior is quite conservative:
        //1. all long nodes ahead should have assignment
        //2. only one should have correct assignment
        //3. this one should have unambiguous path backward to the vertex maybe stopping one link away
        let long_ahead: Vec<Vertex> = self.bounded_dfs(v);

        //println!("Long ahead: {}", long_ahead.iter().map(|x| self.g.v_str(*x)).collect::<Vec<String>>().join(";"));

        if long_ahead.iter().all(|x| self.assignments.is_definite(x.node_id)) {
            let potential_ext: Vec<Vertex> = long_ahead.into_iter()
                .filter(|x| self.assignments.get(x.node_id).unwrap().group == group)
                .collect();
            if potential_ext.len() == 1 {
                let potential_ext = *potential_ext.first().unwrap();
                let mut p = self.unambig_path_forward(potential_ext.rc(), group);
                if !p.in_path(v.node_id) {
                    p = self.try_link(p, v.rc());
                }
                if !p.in_path(v.node_id) {
                    p = self.try_link_with_vertex(p, v.rc(), group);
                }
                if p.trim_to(&v.rc()) > 0 {
                    assert!(p.len() > 1);
                    let p = p.reverse_complement();
                    //println!("Successful jump, path {}", p.print(self.g));
                    return Some(p);
                }
            }
        }

        //println!("Can't jump");

        None
    }

    fn unambig_path_forward(&self, v: Vertex, group: TrioGroup) -> HaploPath {
        let mut p = HaploPath::new(v);
        self.unambig_grow_forward(&mut p, group);
        p
    }

    fn unambig_grow_forward(&self, path: &mut HaploPath, group: TrioGroup) -> usize {
        let mut v = *path.end().unwrap();
        let mut steps = 0;
        while let Some(l) = self.unambiguous_extension(v) {
            let w = l.end;
            if path.in_path(w.node_id) || self.incompatible_assignment(w.node_id, group) {
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
        if let Some(assign) = self.assignments.get(node_id) {
            if TrioGroup::incompatible(assign.group, target_group) {
                return true;
            }
        }
        false
    }

    //maybe move to graph or some GraphAlgoHelper?
    fn unambiguous_extension(&self, v: Vertex) -> Option<Link> {
        //TODO simplify?
        match self.g.outgoing_edge_cnt(v) {
            1 => Some(*self.g.outgoing_edges(v).first().unwrap()),
            _ => None,
        }
    }

}

pub struct HaploPath {
    v_storage: Vec<Vertex>,
    l_storage: Vec<Link>,
    initial_node: usize,
}

//FIXME rename, doesn't know about haplotype!
impl HaploPath {

    fn new(init_v: Vertex) -> HaploPath {
        HaploPath {
            v_storage: vec![init_v],
            l_storage: Vec::new(),
            initial_node: init_v.node_id,
        }
    }

    fn vertices(&self) -> &Vec<Vertex> {
        &self.v_storage
    }

    fn start(&self) -> Option<&Vertex> {
        self.v_storage.first()
    }

    fn end(&self) -> Option<&Vertex> {
        self.v_storage.last()
    }

    fn len(&self) -> usize {
        self.v_storage.len()
    }

    fn links(&self) -> &Vec<Link> {
        &self.l_storage
    }

    fn reverse_complement(self) -> HaploPath {
        //TODO optimize since consuming self
        HaploPath {
            v_storage: self.v_storage.iter().rev().map(|v| v.rc()).collect(),
            l_storage: self.l_storage.iter().rev().map(|l| l.rc()).collect(),
            initial_node: self.initial_node,
        }
    }

    fn trim_to(&mut self, v: &Vertex) -> usize {
        while self.v_storage.len() > 0 && &self.v_storage[self.v_storage.len() - 1] != v {
            self.v_storage.pop();
            self.l_storage.pop();
        }
        self.len()
    }

    fn append(&mut self, l: Link) {
        assert!(self.v_storage.last().unwrap() == &l.start);
        //TODO disable expensive assert?
        assert!(!self.in_path(l.end.node_id));
        self.v_storage.push(l.end);
        self.l_storage.push(l);
    }

    fn in_path(&self, node_id: usize) -> bool {
        self.v_storage.iter().any(|v| v.node_id == node_id)
    }

    fn can_merge_in(&self, path: &HaploPath) -> bool {
        assert!(self.v_storage.last().unwrap() == path.v_storage.first().unwrap());
        !path.v_storage.iter().skip(1).any(|v| self.in_path(v.node_id))
    }

    fn merge_in(&mut self, path: HaploPath) {
        assert!(self.can_merge_in(&path));
        for l in path.l_storage {
            self.append(l);
        }
    }

    pub fn print(&self, g: &Graph) -> String {
        self.v_storage.iter().map(|&v| g.v_str(v)).collect::<Vec<String>>().join(",")
    }

}
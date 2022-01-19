use crate::graph::*;
use std::cmp;
use std::collections::HashSet;
use std::collections::HashMap;
use log::debug;

type DistRange = (usize, usize);

fn shift_range((min, max): DistRange, s: usize) -> DistRange {
    (min + s, max + s)
}

fn merge_range((min1, max1): DistRange, (min2, max2): DistRange) -> DistRange {
    (cmp::min(min1, min2), cmp::max(max1, max2))
}

pub struct Superbubble<'a> {
    g: &'a Graph,
    start_vertex: Vertex,
    end_vertex: Option<Vertex>,
    //vertex to path length range
    reached_vertices: HashMap<Vertex, DistRange>,
}

impl <'a> Superbubble<'a> {

    fn link_dist_range(&self, l: Link) -> DistRange {
        let &r = self.reached_vertices.get(&l.start).unwrap();
        let enode_len = self.g.node(l.end.node_id).length;
        assert!(enode_len >= l.overlap);
        shift_range(r, enode_len - l.overlap)
    }

    pub fn longest_path(&self) -> Path {
        let mut v = self.end_vertex.unwrap();
        let mut longest_dist = self.reached_vertices.get(&v).unwrap().1;
        let mut rc_p = Path::new(v.rc());
        'outer: while v != self.start_vertex {
            //let l = self.heaviest_backtrace.get(v).unwrap();
            for l in self.g.incoming_edges(v) {
                if self.link_dist_range(l).1 == longest_dist {
                    assert!(l.end == v);
                    rc_p.append(l.rc());
                    v = l.start;
                    longest_dist = self.reached_vertices.get(&l.start).unwrap().1;
                    continue 'outer;
                }
            }
            assert!(false);
        }
        rc_p.reverse_complement()
    }

    pub fn shortest_path(&self) -> Path {
        let mut v = self.end_vertex.unwrap();
        let mut shortest_dist = self.reached_vertices.get(&v).unwrap().0;
        let mut rc_p = Path::new(v.rc());
        'outer: while v != self.start_vertex {
            //let l = self.heaviest_backtrace.get(v).unwrap();
            for l in self.g.incoming_edges(v) {
                if self.link_dist_range(l).0 == shortest_dist {
                    assert!(l.end == v);
                    rc_p.append(l.rc());
                    v = l.start;
                    shortest_dist = self.reached_vertices.get(&l.start).unwrap().0;
                    continue 'outer;
                }
            }
            assert!(false);
        }
        rc_p.reverse_complement()
    }

    pub fn vertices(&self) -> impl Iterator<Item=&Vertex> + '_ {
        self.reached_vertices.keys()
    }

    pub fn start_vertex(&self) -> Vertex {
        self.start_vertex
    }

    pub fn end_vertex(&self) -> Vertex {
        self.end_vertex.unwrap()
    }

    //Maybe return Option?
    pub fn length_range(&self) -> (usize, usize) {
        match self.end_vertex {
            None => (0, 0),
            Some(v) => *self.reached_vertices.get(&v).unwrap(),
        }
    }
}

//TODO can be heavily optimized (e.g. no maps, sets, etc)
//TODO support other weights -- currently using max length
//Maybe update to pseudo-code from miniasm paper?
pub struct SbSearchParams {
    pub max_length: usize,
    pub max_diff: usize,
    pub max_count: usize,
}

impl SbSearchParams {
    //all usize values should probably default to max values
    //FIXME provide builder
    pub fn unrestricted() -> SbSearchParams {
        SbSearchParams {
            max_length: usize::MAX,
            max_diff: usize::MAX,
            max_count: usize::MAX,
        }
    }
}

//TODO handle case when first/last vertex have other outgoing/incoming edges
//last vertex case is almost handled
//returns true if no thresholds exceeded
pub fn find_superbubble<'a>(g: &'a Graph, v: Vertex, params: &SbSearchParams) -> Option<Superbubble<'a>> {
    let mut bubble = Superbubble {
        g,
        start_vertex: v,
        reached_vertices: HashMap::new(),
        end_vertex: None,
    };
    if g.outgoing_edge_cnt(bubble.start_vertex) < 2
        //same check, but excluding loops
        || g.outgoing_edges(bubble.start_vertex).iter().filter(|l| l.start != l.end).count() < 2 {
        return None;
    }

    debug!("Adding starting vertex {} to stack", g.v_str(bubble.start_vertex));
    //vertices with all incoming edges considered (can be processed)
    let mut can_be_processed: Vec<Vertex> = vec![bubble.start_vertex];
    bubble.reached_vertices.insert(bubble.start_vertex, (0, 0));

    //reached vertices that can't be processed yet
    let mut not_ready_cnt = 0;
    let mut remaining_incoming: HashMap<Vertex, usize> = HashMap::new();

    while !can_be_processed.is_empty() {
        if bubble.reached_vertices.len() > params.max_count {
            return None;
        }

        let v = can_be_processed.pop().unwrap();
        debug!("Adding vertex {} to the bubble", g.v_str(v));

        if g.outgoing_edge_cnt(v) == 0 {
            debug!("Hit dead-end");
            return None;
        }

        debug!("Looking at neighbors");
        for l in g.outgoing_edges(v) {
            let w = l.end;
            if w == bubble.start_vertex {
                return None;
                //FIXME re-enable after dealing with usage wrt start/end symmetry absense
                //if v != self.start_vertex {
                //    //no loops involiving the start vertex
                //    return false;
                //} else {
                //    //unless self-loop
                //    continue;
                //}
            }

            if !bubble.reached_vertices.contains_key(&w) {
                if bubble.reached_vertices.contains_key(&w.rc()) {
                    debug!("Reverse-complement vertex {} was already reached",
                        g.v_str(w.rc()));
                    return None;
                }
                not_ready_cnt += 1;
                remaining_incoming.insert(w, g.incoming_edge_cnt(w));
                bubble.reached_vertices.insert(w, bubble.link_dist_range(l));
            }
            let rem_inc = remaining_incoming.get_mut(&w).unwrap();
            *rem_inc -= 1;
            //self.reached_vertices.get(w) =
            bubble.reached_vertices.insert(w,
                merge_range(*bubble.reached_vertices.get(&w).unwrap(), bubble.link_dist_range(l)));

            if *remaining_incoming.get(&w).unwrap() == 0 {
                can_be_processed.push(w);
                not_ready_cnt -= 1;
            }
        }

        if can_be_processed.len() == 1 && not_ready_cnt == 0 {
            //FIXME second case is not a classic one, check that it works!
            //Also needs more work to get final vertex!!!
            //|| (can_be_processed.len() == 0 && not_ready_cnt == 1)
            //process last vertex?
            debug!("End node found! Vertex {}", g.v_str(v));
            let v = can_be_processed.pop().unwrap();

            let &(min_len, max_len) = bubble.reached_vertices.get(&v).unwrap();

            let v_len = g.node(v.node_id).length;

            //FIXME it seems like only start_pos is ever checked
            if min_len > v_len && (min_len - v_len) > params.max_length {
                debug!("Length of minimal additional sequence {} exceeded limit {}",
                    min_len - v_len, params.max_length);
                return None;
            }
            if max_len - min_len > params.max_diff {
                debug!("Minimal and maximal lengths differed by {} exceeded limit {}",
                    max_len - min_len, params.max_diff);
                return None;
            }
            bubble.end_vertex = Some(v);
            return Some(bubble);
        }
    }

    debug!("No more nodes could be added");
    debug!("Finished search for starting vertex {}", g.v_str(bubble.start_vertex));
    None
}

type BubbleChain<'a> = Vec<Superbubble<'a>>;

pub fn find_chain_ahead<'a>(g: &'a Graph, init_v: Vertex, params: &SbSearchParams) -> BubbleChain<'a> {
    let mut chain = Vec::new();
    //FIXME no need to check here, since we are marking everything, but useful for general code
    let mut v = init_v;

    loop {
        match find_superbubble(g, v, &params) {
            None => break,
            Some(bubble) => {
                v = bubble.end_vertex();
                chain.push(bubble);
                if v == init_v {
                    break;
                }
            }
        }
    }
    chain
}
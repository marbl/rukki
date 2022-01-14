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

//TODO can be heavily optimized (e.g. no maps, sets, etc)
//TODO support other weights -- currently using max length
//Maybe update to pseudo-code from miniasm paper?
pub struct SuperbubbleFinder<'a> {
    g: &'a Graph,
    start_vertex: Vertex,
    max_length: usize,
    max_diff: usize,
    max_count: usize,

    //FIXME maybe create superbubble struct?
    cnt_: usize,
    //TODO think of alternative definitions of weight (currently: total k-mer multiplicity)
    //vertex to heaviest path weight / path length range
    reached_vertices: HashMap<Vertex, DistRange>,
    end_vertex: Option<Vertex>,
}

impl<'a> SuperbubbleFinder<'a> {

    fn no_link_to_start(&self, v: Vertex) -> bool {
        !self.g.outgoing_edges(v).iter().any(|l| l.end == self.start_vertex)
    }

    fn link_dist_range(&self, l: Link) -> DistRange {
        let &r = self.reached_vertices.get(&l.start).unwrap();
        let enode_len = self.g.node(l.end.node_id).length;
        assert!(enode_len >= l.overlap);
        shift_range(r, enode_len - l.overlap)
    }

    //all usize values should probably default to max values
    //FIXME what is the conventional way to set/change defaults?
    pub fn new(g: &'a Graph, v: Vertex, max_length: usize, max_diff: usize, max_count: usize) -> SuperbubbleFinder {
        SuperbubbleFinder {
            g,
            start_vertex: v,
            max_length,
            max_diff,
            max_count,
            cnt_: 0,
            reached_vertices: HashMap::new(),
            end_vertex: None,
        }
    }

    //TODO handle case when first/last vertex have other outgoing/incoming edges
    //last vertex case is almost handled
    //returns true if no thresholds exceeded
    pub fn find_superbubble(&mut self) -> bool {
        if self.g.outgoing_edge_cnt(self.start_vertex) < 2
            //same check, but excluding loops
            || self.g.outgoing_edges(self.start_vertex).iter().filter(|l| l.start != l.end).count() < 2 {
            return false;
        }

        debug!("Adding starting vertex {} to stack", self.g.v_str(self.start_vertex));
        //vertices with all incoming edges considered (can be processed)
        let mut can_be_processed: Vec<Vertex> = vec![self.start_vertex];
        self.reached_vertices.insert(self.start_vertex, (0, 0));

        //reached vertices that can't be processed yet
        let mut not_ready_cnt = 0;
        //let update_can_be_processed = |v| {
        //    for l in self.g.outgoing_links(v) {
        //        let w = l.end;
        //        //debug!("Considering neighbor {}", self.g.v_str(w));
        //        //handling self-loop at start vertex
        //        if w == self.start_vertex {
        //            //for all other nodes it can't happen at this point
        //            assert!(v == self.start_vertex);
        //            continue;
        //        }
        //        assert!(!self.reached_vertices.contains(neighbour_v));
        //        //debug!("Adding vertex {} to border", self.g.v_str(neighbour_v));
        //        border.insert(neighbour_v);
        //        if (self.check_can_be_processed(neighbour_v)) {
        //            debug!("Adding vertex {} to 'can be processed' set", self.g.v_str(neighbour_v));
        //            can_be_processed.insert(neighbour_v);
        //        }
        //    }
        //};
        //update_can_be_processed(start_vertex_, can_be_processed, border);

        let mut remaining_incoming: HashMap<Vertex, usize> = HashMap::new();

        while !can_be_processed.is_empty() {
            if self.reached_vertices.len() > self.max_count {
                return false;
            }

            let v = can_be_processed.pop().unwrap();
            debug!("Adding vertex {} to the bubble", self.g.v_str(v));

            if self.g.outgoing_edge_cnt(v) == 0 {
                debug!("Hit dead-end");
                return false;
            }

            debug!("Looking at neighbors");
            for l in self.g.outgoing_edges(v) {
                let w = l.end;
                if w == self.start_vertex {
                    if v != self.start_vertex {
                        //no loops involiving the start vertex
                        return false;
                    } else {
                        //unless self-loop
                        continue;
                    }
                }

                if !self.reached_vertices.contains_key(&w) {
                    if self.reached_vertices.contains_key(&w.rc()) {
                        debug!("Reverse-complement vertex {} was already reached",
                            self.g.v_str(w.rc()));
                        return false;
                    }
                    not_ready_cnt += 1;
                    remaining_incoming.insert(w, self.g.incoming_edge_cnt(w));
                    self.reached_vertices.insert(w, self.link_dist_range(l));
                }
                let rem_inc = remaining_incoming.get_mut(&w).unwrap();
                *rem_inc -= 1;
                //self.reached_vertices.get(w) =
                self.reached_vertices.insert(w,
                    merge_range(*self.reached_vertices.get(&w).unwrap(), self.link_dist_range(l)));

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
                debug!("End node found! Vertex {}", self.g.v_str(v));
                let v = can_be_processed.pop().unwrap();

                let &(min_len, max_len) = self.reached_vertices.get(&v).unwrap();

                let v_len = self.g.node(v.node_id).length;

                //FIXME it seems like only start_pos is ever checked
                if min_len > v_len && (min_len - v_len) > self.max_length {
                    debug!("Length of minimal additional sequence {} exceeded limit {}",
                        min_len - v_len, self.max_length);
                    return false;
                }
                if max_len - min_len > self.max_diff {
                    debug!("Minimal and maximal lengths differed by {} exceeded limit {}",
                        max_len - min_len, self.max_diff);
                    return false;
                }
                self.end_vertex = Some(v);
                return true;
            }
        }

        debug!("No more nodes could be added");
        debug!("Finished search for starting vertex {}", self.g.v_str(self.start_vertex));
        false
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

    pub fn end_vertex(&self) -> Option<Vertex> {
        self.end_vertex
    }

    //Maybe return Option?
    pub fn length_range(&self) -> (usize, usize) {
        match self.end_vertex {
            None => (0, 0),
            Some(v) => *self.reached_vertices.get(&v).unwrap(),
        }
    }

}
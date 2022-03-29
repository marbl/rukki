//use rukki::*;
//use rukki::graph_algos::scc;
//use std::fs;
//use std::fs::File;
//use std::io::Write;
////FIXME populate with small corner cases.

//#[test]
//fn manual_tmp_test() {
//    let in_file = "";
//    let out_file = "";
//    let scc_out_file = "";
//    let g = Graph::read(&fs::read_to_string(in_file).unwrap());
//    let sccs = scc::strongly_connected(&g);
//    let (cond, _v_map) = scc::condensation(&g, &sccs, false);
//    let mut output = File::create(out_file).unwrap();
//    write!(output, "{}", cond.as_gfa()).unwrap();
//    let mut output = File::create(scc_out_file).unwrap();
//    for (scc_id, scc) in sccs.iter().enumerate() {
//        write!(output, "scc_{}: {}\n", scc_id, scc.iter().map(|&w| g.v_str(w)).collect::<Vec<String>>().join(",")).unwrap();
//    }
//}
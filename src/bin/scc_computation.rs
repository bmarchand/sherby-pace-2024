//use clap::Parser;
use sherby_pace_2024::*;

fn main() {
//    let args = Cli::parse();

    let graph_name = std::env::args().nth(1).expect("expecting a path to a graph file.");
    let graph: Graph = parse_graph(&graph_name);

    let mut crossing_dict = orientable_crossing_values(&graph);

    // twins pre-processing
    let twin_mapping = find_twins(&graph);
    println!(
        "{:?} twins out of {:?} vertices",
        twin_mapping.len(),
        graph.bnodes.len()
    );

    let graph = process_twins_graph(graph, &twin_mapping);
    println!(
        "num vertices after twin processing {:?}",
        graph.bnodes.len()
    );

    crossing_dict = process_twins_crossing_dict(&twin_mapping, &crossing_dict);
    // end twins pre-processing

    let _sccs = compute_scc(&graph, &crossing_dict);
    //    let mut h = petgraph::graph::Graph::<usize, usize>::new();
    //
    //    let mut map: HashMap<usize, petgraph::graph::NodeIndex> = HashMap::new();
    //
    //    for u in &graph.bnodes {
    //        let nu = h.add_node(u.id);
    //        map.insert(u.id, nu);
    //    }
    //
    //    for (u, v) in crossing_dict.keys() {
    //        if let Some(c) = crossing_dict.get(&(*v, *u)) {
    //            if *c > *crossing_dict.get(&(*u, *v)).unwrap() {
    //                h.add_edge(*map.get(u).unwrap(), *map.get(v).unwrap(), 1);
    //            } else if *c < *crossing_dict.get(&(*u, *v)).unwrap() {
    //                h.add_edge(*map.get(v).unwrap(), *map.get(u).unwrap(), 1);
    //            }
    //        }
    //    }
    //
    //    for u in &graph.bnodes {
    //        for v in &graph.bnodes {
    //            if u.id != v.id {
    //                if u.right <= v.left {
    //                    h.add_edge(*map.get(&u.id).unwrap(), *map.get(&v.id).unwrap(), 1);
    //                }
    //            }
    //        }
    //    }
    //
    //    let scc = tarjan_scc(&h);
    //    println!("number of independent components {:?}", scc.len());
    //    for c in &scc {
    //        println!("{:?}", c);
    //    }
}

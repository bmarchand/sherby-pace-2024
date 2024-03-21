use clap::Parser;
use peak_alloc::PeakAlloc;
use sherby_pace_2024::*;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    let args = Cli::parse();

    let graph: Graph = parse_graph(&args.graph);

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

    // scc computations
    let sccs = compute_scc(&graph, &crossing_dict);
    println!("number of SCCs {:?}", sccs.len());

    // final solution init
    let mut vec: Vec<usize> = Vec::new();

    //    for scc in &sccs {
    //        print_edges(scc);
    //    }

    // main calls
    for scc in &sccs {
        let vec_scc = kobayashi_tamaki(&scc, &crossing_dict).unwrap();
        vec.extend_from_slice(&vec_scc);
    }

    // twins post-processing
    let vec = add_twins(vec, &twin_mapping);
    // end twin post-processing

    // Writing result in output file (name same as input, extension changed)
    let mut outname = args.graph.clone();
    outname.set_extension("sol");
    let v: Vec<String> = vec.into_iter().map(|x| x.to_string()).collect();
    let _ = std::fs::write(outname, v.join("\n"));

    let peak_mem = PEAK_ALLOC.peak_usage_as_mb();
    println!("peak memory: {} mb", peak_mem);
}

//fn mergeable_adjacent_fixed_vertices(graph: &Graph) -> Vec<(usize,usize)> {
//
//    let mergeable_pairs: Vec<(usize,usize)> = Vec::new();
//
//    // #1.5 prep work for orientable pairs
//    let mut left_end: HashMap<usize, usize> = HashMap::new();
//    let mut right_end: HashMap<usize, usize> = HashMap::new();
//    let mut neighbors: HashMap<usize, Vec<usize>> = HashMap::new();
//    for u in &graph.bnodes {
//        left_end.insert(u.id, u.left);
//        right_end.insert(u.id, u.right);
//        neighbors.insert(u.id, u.neighbors.clone());
//    }
//
//    let mut orientable_pairs: HashSet<(usize, usize)> = HashSet::new();
//    let mut active_vertices: HashSet<usize> = HashSet::new();
//    // #2 computing orientable pairs
//    for k in 0..(&graph.anodes.len()-1) {
//        let x = graph.anodes[k];
//        let y = graph.anodes[k+1];
//        for u in &x.neighbors {
//            if x.id == *left_end.get(u).unwrap() {
//                active_vertices.insert(*u);
//            }
//            if x.id == *right_end.get(u).unwrap() {
//                active_vertices.remove(u);
//            }
//        }
//        let mut mergeable = true;
//
//        // no vertex active at both should be connected to both
//        for u in active_vertices {
//            if u.neighbors.contains(x) && u.neighbors.contains(y) {
//                mergeable = false;
//                break
//            }
//        }
//
//        if mergeable {
//            for u in active_vertices {
//                for v in active_vertices {
//                    if u.neighbors.contains(x) && v.neighbors.contains(y) {
//                        mergeable = false;
//                        break;
//                    }
//                }
//            }
//        }
//
//        if mergeable {
//            mergeable_pairs.push((,))
//        }
//
//        for u in &y.neighbors {
//            if y.id == *right_end.get(u).unwrap() {
//                active_vertices.remove(u);
//            }
//        }
//
//    }
//    return mergeable_pairs;
//}

//fn print_edges(graph: &Graph) {
//
//    let mut num_edges = 0;
//    for u in &graph.bnodes {
//        for v in &u.neighbors {
//            num_edges += 1;
//        }
//    }
//
//    println!("p ocr ? {:?} {:?}", graph.bnodes.len(), num_edges);
//
//    for u in &graph.bnodes {
//        for v in &u.neighbors {
//            println!("{:?} {:?}",v,u.id);
//        }
//    }
//}

//fn process_twins_ints(twin_mapping: &HashMap<usize,usize>,ints: Vec<(usize, usize, usize)>) -> Vec<(usize,usize,usize)> {
//
//    let mut new_ints: Vec<(usize,usize,usize)> = Vec::new();
//
//    for tup in &ints {
//        if twin_mapping.contains_key(&tup.0) {
//            continue
//        }
//        new_ints.push(*tup);
//    }
//    return new_ints;
//}

///// Computing the number of binary numbers having Hamming
///// weight w and L bits in their decomposition (with potential
///// 0s on the left side)
//fn num_fixed_hamming_weigth(L: usize, w: usize) -> u64 {
//
//}
//
//fn subset_rank(subset: Vec<bool>) -> u64 {
//
//
//}

///// Computation of all crossing values c(u,v) for u,v two
///// vertices of Y, the side of the graph that must be ordered.
///// The result is an Hashmap associating keys (u,v) to
///// the integer value c(u,v), the number of crossings
///// obtained between edges adjacent to u and edges
///// adjacent to v if u is placed before v in the order.
/////
///// It works by first compuring $d^{<x}(u)$, the number of neighbors
///// of $u$ (a vertex of B) that are strictly before x, a vertex in A.
///// Then, it uses the formula $$c(u,v) = \sum_{x\in N(u)} d^{<x}(v)$$
//fn crossing_values(graph: &Graph) -> HashMap<(usize, usize), usize> {
//
//    let mut d_less_than_x: HashMap<(usize,usize), usize> = HashMap::new();
//
//    // #1 compute all d^<x(u) O(|X||Y|)
//    for u in &graph.bnodes {
//        let mut d: usize = 0;
//        for x in &graph.anodes {
//            d_less_than_x.insert((u.id,x.id),d);
//            if u.neighbors.contains(&x.id) {
//                d += 1;
//            }
//        }
//    }
//
//    // #2 compute all crossing values
//    let mut crossing_dict: HashMap<(usize,usize), usize> = HashMap::new();
//
//    for u in &graph.bnodes {
//        for v in &graph.bnodes {
//            if u.id==v.id {
//                continue;
//            }
//            for x in &u.neighbors {
//                if let Some(c) = crossing_dict.get(&(u.id,v.id)) {
//                    crossing_dict.insert((u.id,v.id),*c+d_less_than_x.get(&(v.id,*x)).unwrap());
//                }
//                else {
//                    crossing_dict.insert((u.id,v.id),*d_less_than_x.get(&(v.id,*x)).unwrap());
//                }
//            }
//        }
//    }
//
//    return crossing_dict;
//}

use clap::Parser;
use peak_alloc::PeakAlloc;
use sherby_pace_2024::*;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    let args = Cli::parse();

    let graph: Graph = parse_graph_cutwidth(&args.graph);

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

    // main calls
    for scc in &sccs {
//        let vec_scc = kobayashi_tamaki(&scc, &crossing_dict).unwrap();
        let vec_scc = recursive_kt(&scc, &crossing_dict).unwrap();
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

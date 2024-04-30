/// Solver for the cutwidth track.
/// Its only difference with the exact track
/// solver (main.rs) is the parsing of
/// the graph, which ignores the ordering
/// that is given.
//use clap::Parser;
//use peak_alloc::PeakAlloc;
use sherby_pace_2024::*;
//use log::info;

//#[global_allocator]
//static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
//    let args = Cli::parse();
    let graph_name = std::env::args().nth(1).expect("expecting a path to a graph file.");
    let solution_name = std::env::args().nth(2).expect("expecting a path to a graph file.");

//    if args.dfas{
//    	//info!("dfas mode activated");
//    }

    let graph: Graph = parse_graph_cutwidth(&graph_name);

    let mut crossing_dict = orientable_crossing_values(&graph);

    // twins pre-processing
    let twin_mapping = find_twins(&graph);
    //info!(
//        "{:?} twins out of {:?} vertices",
//        twin_mapping.len(),
//        graph.bnodes.len()
//    );

    let graph = process_twins_graph(graph, &twin_mapping);
    //info!(
//        "num vertices after twin processing {:?}",
//        graph.bnodes.len()
//    );

    crossing_dict = process_twins_crossing_dict(&twin_mapping, &crossing_dict);
    // end twins pre-processing

    // scc computations
    let sccs = compute_scc(&graph, &crossing_dict);
    //info!("number of SCCs {:?}", sccs.len());

    // final solution init
    let mut vec: Vec<usize> = Vec::new();

    // main calls
	let mut cptscc : usize = 1;
    for scc in &sccs {
//        let vec_scc = kobayashi_tamaki(&scc, &crossing_dict).unwrap();
        //
                let instance_size = total_instance_size(&scc);
        
		if instance_size > 10000000 && scc.bnodes.len() > 1{
			//info!("Solving scc #{}, size={}", cptscc, scc.bnodes.len());
			let (_cost, vec_scc) = solve_dfas( &scc, &crossing_dict, cptscc );
			vec.extend_from_slice(&vec_scc);
		}
		else{
			//let vec_scc = recursive_kt(&scc, &crossing_dict).unwrap();
			let vec_scc = kobayashi_tamaki(&scc, &crossing_dict).unwrap();
			vec.extend_from_slice(&vec_scc);
		}
		
		cptscc += 1;
    }

    // twins post-processing
    let vec = add_twins(vec, &twin_mapping);
    // end twin post-processing

    // Writing result in output file (name same as input, extension changed)
    let outname = solution_name.clone();
//    outname.set_extension("sol");
    let v: Vec<String> = vec.into_iter().map(|x| x.to_string()).collect();
    let _ = std::fs::write(outname, v.join("\n"));

//    let peak_mem = PEAK_ALLOC.peak_usage_as_mb();
    //info!("peak memory: {} mb", peak_mem);
}

use ::std::collections::HashMap;
use clap::Parser;
use peak_alloc::PeakAlloc;
use sherby_pace_2024::*;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    let args = Cli::parse();

    let graph: Graph = parse_graph(&args.graph);

    let crossing_dict = crossing_values(&graph);

    let vec = kobayashi_tamaki(&graph, &crossing_dict).unwrap();

    // Writing result in output file (name same as input, extension changed)
    let mut outname = args.graph.clone();
    outname.set_extension("no_preprocess_sol");
    let v: Vec<String> = vec.into_iter().map(|x| x.to_string()).collect();
    let _ = std::fs::write(outname, v.join("\n"));

    let peak_mem = PEAK_ALLOC.peak_usage_as_mb();
    println!("peak memory: {} mb", peak_mem);
}

/// Computation of all crossing values c(u,v) for u,v two
/// vertices of Y, the side of the graph that must be ordered.
/// The result is an Hashmap associating keys (u,v) to
/// the integer value c(u,v), the number of crossings
/// obtained between edges adjacent to u and edges
/// adjacent to v if u is placed before v in the order.
///
/// It works by first compuring $d^{<x}(u)$, the number of neighbors
/// of $u$ (a vertex of B) that are strictly before x, a vertex in A.
/// Then, it uses the formula $$c(u,v) = \sum_{x\in N(u)} d^{<x}(v)$$
fn crossing_values(graph: &Graph) -> HashMap<(usize, usize), usize> {
    let mut d_less_than_x: HashMap<(usize, usize), usize> = HashMap::new();

    // #1 compute all d^<x(u) O(|X||Y|)
    for u in &graph.bnodes {
        let mut d: usize = 0;
        for x in &graph.anodes {
            d_less_than_x.insert((u.id, x.id), d);
            if u.neighbors.contains(&x.id) {
                d += 1;
            }
        }
    }

    // #2 compute all crossing values
    let mut crossing_dict: HashMap<(usize, usize), usize> = HashMap::new();

    for u in &graph.bnodes {
        for v in &graph.bnodes {
            if u.id == v.id {
                continue;
            }
            for x in &u.neighbors {
                if let Some(c) = crossing_dict.get(&(u.id, v.id)) {
                    crossing_dict
                        .insert((u.id, v.id), *c + d_less_than_x.get(&(v.id, *x)).unwrap());
                } else {
                    crossing_dict.insert((u.id, v.id), *d_less_than_x.get(&(v.id, *x)).unwrap());
                }
            }
        }
    }

    return crossing_dict;
}

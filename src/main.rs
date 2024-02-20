use clap::Parser;
use std::collections::HashMap;

#[derive(Parser)]
struct Cli {
    /// The path to the graph
    graph: std::path::PathBuf,
}

#[derive(Default, Debug)]
struct BNode {
    id: usize,
    /// List of neighbors
    neighbors: Vec<usize>,
    /// Left-most neighbor
    left: usize,
    /// Right-most neighbor
    right: usize,
}

#[derive(Default, Debug)]
struct Graph {
    /// vector of BNodes
    bnodes: Vec<BNode>,
}

fn main() {
    let args = Cli::parse();

    let content = std::fs::read_to_string(&args.graph).expect("could not read file");

    // graph initialization
    let mut graph: Graph = Default::default();
    let mut n0: usize = 0;

    // parsing file
    for line in content.lines() {
        if line.starts_with("p") {
            n0 = line
                .split_whitespace()
                .nth(2)
                .expect("should be the number of vertices in A")
                .parse()
                .unwrap();
            let n1: usize = line
                .split_whitespace()
                .nth(3)
                .expect("should be the number of vertices in B")
                .parse()
                .unwrap();
            for y in (n0 + 1)..=(n0 + n1) {
                graph.bnodes.push(BNode {
                    id: y,
                    left: n0,
                    right: 0,
                    ..Default::default()
                });
            }
            continue;
        }
        println!("{}", line);

        let a: usize = line
            .split_whitespace()
            .nth(0)
            .expect("should be an A vertex")
            .parse()
            .unwrap();

        let b: usize = line
            .split_whitespace()
            .nth(1)
            .expect("should be a B vertex")
            .parse()
            .unwrap();

        graph.bnodes[b - n0 - 1].neighbors.push(a);
        if graph.bnodes[b - n0 - 1].left > a {
            graph.bnodes[b - n0 - 1].left = a;
        }
        if graph.bnodes[b - n0 - 1].right < a {
            graph.bnodes[b - n0 - 1].right = a;
        }
    }
    println!("graph {:?}", graph);

    let ints = nice_interval_repr(&graph);
    println!("nice interval representation {:?}", ints);

    // TODO: implement Kobayashi-Tamaki algo
    let crossing_dict = crossing_values(&graph);

    println!("crossing values {:?}", crossing_dict);

    let mut k: u32 = 0;
    let mut outname = args.graph.clone();
    outname.set_extension("sol");

    loop {
        let res = kobayashi_tamaki(&ints, &crossing_dict, k);

        println!("current k value: {k}");

        match res {
            Ok(vec) => {
                let v: Vec<String> = vec.into_iter().map(|x| x.to_string()).collect();

                let _ = std::fs::write(outname, v.join("\n"));
                break;
            }
            Err(_) => {
                k += 1;
                continue;
            }
        }
    }
}

fn kobayashi_tamaki(
    ints: &Vec<(usize, usize, usize)>,
    crossing_dict: &HashMap<(usize, usize), u32>,
    k: u32,
) -> Result<Vec<usize>, String> {
    // This is a dummy implem, whose only purpose is to match final behavior
    if k < 2 {
        return Err("k not high enough".to_string());
    }

    let mut v = Vec::new();
    for p in ints {
        v.push(p.0);
    }
    Ok(v)
}

fn crossing_values(graph: &Graph) -> HashMap<(usize, usize), u32> {
    // This implementation is simple but NOT OPTIMAL.

    let mut crossing_dict = HashMap::new();

    for node1 in &graph.bnodes {
        for node2 in &graph.bnodes {
            let mut c: u32 = 0;
            if node1.id != node2.id {
                for x1 in &node1.neighbors {
                    for x2 in &node2.neighbors {
                        if node1.id < node2.id && x1 > x2 || node1.id > node2.id && x2 < x1 {
                            c += 1;
                        }
                    }
                }
                crossing_dict.insert((node1.id, node2.id), c);
            }
        }
    }

    return crossing_dict;
}

fn exclude_first(p: &(usize, i32, i32)) -> (i32, i32) {
    (p.1, p.2)
}

fn nice_interval_repr(graph: &Graph) -> Vec<(usize, usize, usize)> {
    let mut p = Vec::new();

    for node in &graph.bnodes {
        let degree: i32 = node.neighbors.len() as i32;
        let pl: (usize, i32, i32) = (node.id, node.left as i32, -degree);
        let pr: (usize, i32, i32) = (node.id, node.right as i32, 2 * degree);
        p.push(pl);
        p.push(pr);
    }

    // lexico-graphic order while ignoring first
    p.sort_by(|a, b| exclude_first(a).cmp(&exclude_first(b)));

    let mut left_end = HashMap::new();
    let mut right_end = HashMap::new();

    for (idx, tup) in p.iter().enumerate() {
        if tup.2 < 0 {
            left_end.insert(tup.0, idx);
        }
        if tup.2 > 0 {
            right_end.insert(tup.0, idx);
        }
    }

    let mut ret = Vec::new();

    for node in &graph.bnodes {
        ret.push((
            node.id,
            *left_end.get(&node.id).unwrap(),
            *right_end.get(&node.id).unwrap(),
        ));
    }

    return ret;
}

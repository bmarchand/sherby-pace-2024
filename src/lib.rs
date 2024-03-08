use clap::Parser;
use std::path::PathBuf;

#[derive(Default, Debug)]
pub struct Graph {
    /// vector of BNodes
    pub bnodes: Vec<BNode>,
    pub anodes: Vec<ANode>,
}

#[derive(Default, Debug)]
pub struct ANode {
    pub id: usize,
    pub neighbors: Vec<usize>,
}

#[derive(Default, Debug)]
pub struct BNode {
    pub id: usize,
    /// List of neighbors
    pub neighbors: Vec<usize>,
    /// Left-most neighbor
    pub left: usize,
    /// Right-most neighbor
    pub right: usize,
}

#[derive(Parser)]
pub struct Cli {
    /// The path to the graph
    pub graph: std::path::PathBuf,
}

pub fn parse_graph(file_name: &PathBuf) -> Graph {
    let content = std::fs::read_to_string(file_name).expect("could not read file");

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
            for x in 1..=n0 {
                graph.anodes.push(ANode {
                    id: x,
                    ..Default::default()
                });
            }
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

        graph.anodes[a-1].neighbors.push(b);
        graph.bnodes[b - n0 - 1].neighbors.push(a);
        if graph.bnodes[b - n0 - 1].left > a {
            graph.bnodes[b - n0 - 1].left = a;
        }
        if graph.bnodes[b - n0 - 1].right < a {
            graph.bnodes[b - n0 - 1].right = a;
        }
    }
    return graph;
}

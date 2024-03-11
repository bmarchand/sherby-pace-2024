use clap::Parser;
use std::collections::HashMap;
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

pub fn process_twins_graph(mut graph: Graph, twin_mapping: &HashMap<usize,usize>) -> Graph {

    let mut new_bnodes: Vec<BNode> = Vec::new();

    for u in graph.bnodes.into_iter() {
        if twin_mapping.contains_key(&u.id) {
            continue
        }
        new_bnodes.push(u);
    }
    
    graph.bnodes = new_bnodes;

    return graph;
}

pub fn find_twins(graph: &Graph) -> HashMap<usize,usize> {
    
    // main hash table
    let mut h: HashMap<Vec<usize>,Vec<usize>> = HashMap::new();

    for u in &graph.bnodes {
        if let Some(v) = h.get_mut(&u.neighbors) {
            (*v).push(u.id);
        }
        else {
            h.insert(u.neighbors.clone(), vec![u.id]);
        }
    }

    // producing mapping
    let mut mapping: HashMap<usize,usize> = HashMap::new();
    
    for (_, v) in h.iter() {
        if v.len() > 1 {
            let m: usize = *v.iter().min().unwrap();
            for u in v {
                if m!=*u {
                    mapping.insert(*u,m);
                }
            }
        }
    }

    return mapping;
}

//pub fn from_twin_list_to_mapping(twin_list: &HashSet<(usize,usize)>) -> HashMap<usize,usize> {
//    // let us first figure out a mapping u->v (u merged into v) for every twin
//    // we want CONSISTENCY: i.e. we want to avoid u merged into v and then w merged into u
//    let mut mapping: HashMap<usize, usize> = HashMap::new();
//
//    // for consistency we build a sorted vector of twins
//    let mut twin_vector: Vec<(usize,usize)> = Vec::new();
//    for (u,v) in twin_list {
//        twin_vector.push((*u,*v));
//    }
//    twin_vector.sort_by(|(u,v),(w,x)| std::cmp::min(w,x).cmp(std::cmp::min(u,v)));
//
//    // final mapping
//    for (u,v) in &twin_vector {
//        let mut w = u;
//        let mut x = v;
//        if let Some(y) = mapping.get(u) {
//            w = y;
//        }
//        if let Some(z) = mapping.get(v) {
//            x = z;
//        }
//        if u < v {
//            mapping.insert(*v,*w);
//        }
//        else {
//            mapping.insert(*u,*x);
//        }
//    }
//    return mapping;
//}

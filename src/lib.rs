use clap::Parser;
use clap_stdin::MaybeStdin;
use petgraph::algo::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::path::PathBuf;
use log::info;

use std::fs::File;
use std::io::Write;
use std::io::{self,BufRead};
use std::process::Command;
use std::path::Path;
use std::fs;


/// The struct (in other languages: class) used
/// to represent an instance (a bipartite graph
/// for which the B side has to ordered while
/// the A side has a fixed order.
#[derive(Default, Debug)]
pub struct Graph {
    /// vector of BNodes
    pub bnodes: Vec<BNode>,
    pub anodes: Vec<ANode>,
}

/// The struct representing a node in the
/// fixed-ordering layer of an instance.
#[derive(Default, Debug)]
pub struct ANode {
    pub id: usize,
    pub neighbors: Vec<usize>,
}

/// The struct representing a node in the
/// free layer of an instance, the one
/// that must be ordered.
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
    pub graph: MaybeStdin<PathBuf>,
    /// The path to the solution
    pub solution: MaybeStdin<PathBuf>,
    
    #[arg(short, long)]
    pub dfas: bool,
    
    #[arg(short, long)]
    pub random: bool,
    
}

/// The same as parse_graph, except that
/// the ordering (i.e. the lines with only one vertex
/// in the cutwidth instances) is ignored.
pub fn parse_graph_cutwidth(file_name: &PathBuf) -> Graph {
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

        if line.split_whitespace().count() < 2 {
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

        graph.anodes[a - 1].neighbors.push(b);
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

/// The fuction parsing a file representing a graph.
/// It simply fills the vectors of A nodes and
/// B nodes of a Graph struct with the required information
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

        graph.anodes[a - 1].neighbors.push(b);
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

/// Given pairs of twin B-vertices (i.e. with the same neighbors),
/// as described by the twin_mapping HashMap, and a graph, returns
/// a new graph in which twins have been merged. A subtlety arises
/// from groups of more than 2 vertices that are all twins.
/// (they are merged into the vertex with the lowest index among the
/// group). On the graph level, merging v into u simply means removing
/// v. The modification of crossing values occurs in its own function.
pub fn process_twins_graph(mut graph: Graph, twin_mapping: &HashMap<usize, usize>) -> Graph {
    let mut new_bnodes: Vec<BNode> = Vec::new();

    for u in graph.bnodes.into_iter() {
        if twin_mapping.contains_key(&u.id) {
            continue;
        }
        new_bnodes.push(u);
    }

    graph.bnodes = new_bnodes;

    return graph;
}

/// Given sets of twins, as represented by
/// the twin_mapping map, merges some entries
/// of the crossing_dict HashMap to account
/// for merges of twins. (one can see
/// the result of merging vertices as
/// the creation of multiedges)
pub fn process_twins_crossing_dict(
    twin_mapping: &HashMap<usize, usize>,
    crossing_dict: &HashMap<(usize, usize), usize>,
) -> HashMap<(usize, usize), usize> {
    let mut new_crossing_dict: HashMap<(usize, usize), usize> = HashMap::new();

    for ((u, v), c) in crossing_dict.iter() {
        let mut w = u;
        let mut x = v;
        if let Some(y) = twin_mapping.get(u) {
            w = y;
        }
        if let Some(z) = twin_mapping.get(v) {
            x = z;
        }
        if let Some(c2) = new_crossing_dict.get(&(*w, *x)) {
            new_crossing_dict.insert((*w, *x), *c2 + *c);
        } else {
            new_crossing_dict.insert((*w, *x), *c);
        }
    }
    return new_crossing_dict;
}

/// Given a graph, returns a HashMap representing sets
/// of twins. In this hashmap, for any two twins
/// u and v such that u < v, there is an entry
/// v -> u. (v will be merged into u)
pub fn find_twins(graph: &Graph) -> HashMap<usize, usize> {
    // below: main hash table,
    // which maps neighborhoods to the list of vertices
    // having this neighborhood.
    let mut h: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();

    for u in &graph.bnodes {
        if let Some(v) = h.get_mut(&u.neighbors) {
            (*v).push(u.id);
        } else {
            h.insert(u.neighbors.clone(), vec![u.id]);
        }
    }

    // producing mapping
    let mut mapping: HashMap<usize, usize> = HashMap::new();

    for (_, v) in h.iter() {
        if v.len() > 1 {
            let m: usize = *v.iter().min().unwrap();
            for u in v {
                if m != *u {
                    mapping.insert(*u, m);
                }
            }
        }
    }

    return mapping;
}

/// Computes crossing values for orientable pairs only.
/// A pair (u,v) is orientable if neither $r_u\leq l_v$ nor $r_v\leq l_u$.
/// It does so by using the same formula as the function crossing_value,
/// except that $d^{<x}(v)$ is only computed for a $x\in [l_v,r_v]$.
/// Then, c(u,v) is only computed for orientable pairs.
/// If a value of $d^{<x}(v)$ is requested for x smaller than $l_v$ or
/// larger than $r_v$ then it is computed on the flies with its trivial values.
pub fn orientable_crossing_values(graph: &Graph) -> HashMap<(usize, usize), usize> {
    let mut d_less_than_x: HashMap<(usize, usize), usize> = HashMap::new();

    // #1 compute d^<x(u) values that matter
    for u in &graph.bnodes {
        let mut d: usize = 0;
        for x in (u.left)..=(u.right) {
            d_less_than_x.insert((u.id, x), d);
            if u.neighbors.contains(&x) {
                d += 1;
            }
        }
    }

    // #1.5 prep work for orientable pairs
    let mut left_end: HashMap<usize, usize> = HashMap::new();
    let mut right_end: HashMap<usize, usize> = HashMap::new();
    let mut neighbors: HashMap<usize, Vec<usize>> = HashMap::new();
    for u in &graph.bnodes {
        left_end.insert(u.id, u.left);
        right_end.insert(u.id, u.right);
        neighbors.insert(u.id, u.neighbors.clone());
    }

    let mut orientable_pairs: HashSet<(usize, usize)> = HashSet::new();
    let mut active_vertices: HashSet<usize> = HashSet::new();
    // #2 computing orientable pairs
    for x in &graph.anodes {
        for u in &x.neighbors {
            if x.id == *left_end.get(u).unwrap() {
                active_vertices.insert(*u);
            }
            if x.id == *right_end.get(u).unwrap() {
                active_vertices.remove(u);
            }
            for w in &active_vertices {
                orientable_pairs.insert((*u, *w));
                orientable_pairs.insert((*w, *u));
            }
        }
    }

    // #3 compute crossing values for orientable pairs
    let mut crossing_dict: HashMap<(usize, usize), usize> = HashMap::new();

    for (u, v) in orientable_pairs {
        if u == v {
            continue;
        }
        for x in neighbors.get(&u).unwrap() {
            // reconstructing d<x(v)
            let mut d = 0;
            if x > right_end.get(&v).unwrap() {
                d = neighbors.get(&v).unwrap().len();
            } else {
                if x > left_end.get(&v).unwrap() {
                    d = *d_less_than_x.get(&(v, *x)).unwrap();
                }
            }

            // adding to crossing value
            if let Some(c) = crossing_dict.get(&(u, v)) {
                crossing_dict.insert((u, v), *c + d);
            } else {
                crossing_dict.insert((u, v), d);
            }
        }
    }

    return crossing_dict;
}

/// Decomposes a graph into its strongly connected components.
/// The directed graph has the free layer B as its vertices,
/// and there is an edge u->v if c(u,v) < c(v,u).
/// Returns a vector of graphs, which can each be treated
/// separately.
pub fn compute_scc(graph: &Graph, crossing_dict: &HashMap<(usize, usize), usize>) -> Vec<Graph> {
    // just mapping id to bnode to remember it
    let mut id_to_bnodes: HashMap<usize, &BNode> = HashMap::new();
    for u in &graph.bnodes {
        id_to_bnodes.insert(u.id, u);
    }

    // directed graph
    let mut h = petgraph::graph::Graph::<usize, usize>::new();

    let mut map: HashMap<usize, petgraph::graph::NodeIndex> = HashMap::new();
    let mut inverse_map: HashMap<usize, usize> = HashMap::new();

    for u in &graph.bnodes {
        let nu = h.add_node(u.id);
        map.insert(u.id, nu);
        inverse_map.insert(nu.index(), u.id);
    }

    for (u, v) in crossing_dict.keys() {
        if let Some(c) = crossing_dict.get(&(*v, *u)) {
            let c2 = crossing_dict.get(&(*u, *v)).unwrap();
            if *c > *c2 {
                let a = *map.get(u).unwrap();
                let b = *map.get(v).unwrap();
                if !h.contains_edge(a, b) {
                    h.add_edge(a, b, *c - *c2);
                }
            } else if *c < *c2 {
                let a = *map.get(v).unwrap();
                let b = *map.get(u).unwrap();
                if !h.contains_edge(a, b) {
                    h.add_edge(a, b, *c2 - *c);
                }
            }
        }
    }

    for u in &graph.bnodes {
        for v in &graph.bnodes {
            if u.id != v.id {
                if u.right <= v.left && u.left < v.right {
                    let a = *map.get(&u.id).unwrap();
                    let b = *map.get(&v.id).unwrap();
                    if !h.contains_edge(a, b) {
                        h.add_edge(a, b, 30000);    //TODO ML: we assume that 30000 is larger than anything
                    }
                }
            }
        }
    }

    info!(
        "h has {:?} edges over {:?} possible",
        h.edge_count(),
        h.node_count() * (h.node_count() - 1) / 2
    );

    let sccs = tarjan_scc(&h);

    let mut graph_vec: Vec<Graph> = Vec::new();
    for scc in &sccs {
        let mut graph_scc: Graph = Default::default();
        for u in scc {
            let id = inverse_map.get(&u.index()).unwrap();
            let bnode = *id_to_bnodes.get(&id).unwrap();
            graph_scc.bnodes.push(BNode {
                id: bnode.id,
                left: bnode.left,
                right: bnode.right,
                neighbors: bnode.neighbors.clone(),
                ..Default::default()
            });
        }
        graph_vec.push(graph_scc);
    }
    graph_vec.reverse();

    return graph_vec;
}

/// Add twins back into a solution computed on the reduced graph.
/// if u was merged into v, then u is reintroduced right next to v
/// in the ordering.
pub fn add_twins(vec: Vec<usize>, twin_mapping: &HashMap<usize, usize>) -> Vec<usize> {
    let mut inverse_mapping: HashMap<usize, Vec<usize>> = HashMap::new();

    for (u, v) in twin_mapping.iter() {
        if let Some(vector) = inverse_mapping.get_mut(v) {
            (*vector).push(*u);
        } else {
            let mut vector = Vec::new();
            vector.push(*u);
            inverse_mapping.insert(*v, vector);
        }
    }

    let mut new_vec: Vec<usize> = Vec::new();
    for u in vec {
        if let Some(vector) = inverse_mapping.get(&u) {
            for v in vector {
                new_vec.push(*v);
            }
        }
        new_vec.push(u);
    }
    return new_vec;
}

/// a simple function that excludes the first component
/// of a tuple, so that it may be ignored in a lexico-graphic
/// ordering
fn exclude_first(p: &(usize, i32, i32)) -> (i32, i32) {
    (p.1, p.2)
}

/// Transforms an integer x into the sub-set
/// of v it represents. if the p-th bit
/// (starting from the right) of x is 1,
/// then the p-th element of v is in the
/// subset.
fn int2vec(x: usize, v: &Vec<usize>) -> Vec<usize> {
    let mut s = Vec::new();
    for (i, y) in v.iter().enumerate() {
        // if i-th bit of x is 1, y in s
        if (x >> i) & 1 != 0 {
            s.push(*y);
        }
    }
    return s;
}

/// Transforms a subset s of a set v into
/// an integer representing it. The p-th
/// bit (starting from the right) of this integer
/// is 1 if the p-th element of v is in s.
fn vec2int(s: &Vec<usize>, v: &Vec<usize>) -> usize {
    let mut m: HashMap<usize, usize> = HashMap::new();
    for (i, y) in v.iter().enumerate() {
        m.insert(*y, i);
    }
    let mut x: usize = 0;
    for y in s {
        x += 1 << m.remove(&y).unwrap();
    }
    return x;
}

/// Computes the numbers of crossings incurred
/// by positioning x after all elements of a set
/// s. i.e. the sum over v in s of c(v,x)
fn set_crossing(s: &Vec<usize>, x: usize, crossing_dict: &HashMap<(usize, usize), usize>) -> usize {
    let mut c: usize = 0;
    for u in s {
        c += crossing_dict.get(&(*u, x)).unwrap_or(&(0 as usize));
    }
    return c;
}

/// Computes a nice interval representation, in the sense of the
/// \mathcal{J} set of intervals in the article of Kobayashi and
/// Tamaki.
///
/// The entries of the output vector are of the form (y, a_y, b_y).
/// I.e. the identifier of the node, the left end of the interval
/// and the right-end of the interval
fn nice_interval_repr(graph: &Graph) -> Vec<(usize, usize, usize)> {
    let mut p = Vec::new();

    for node in &graph.bnodes {
        let degree: i32 = node.neighbors.len() as i32;
        let pl: (usize, i32, i32) = (node.id, node.left as i32, 10 * (degree - 1));
        let pr: (usize, i32, i32) = (node.id, node.right as i32, -10 * (degree - 1) + 1);
        p.push(pl);
        p.push(pr);
    }

    // lexico-graphic order while ignoring first
    p.sort_by(|a, b| exclude_first(a).cmp(&exclude_first(b)));

    let mut left_end = HashMap::new();
    let mut right_end = HashMap::new();

    for (idx, tup) in p.iter().enumerate() {
        if tup.2 >= 10 || tup.2 == 0 {
            left_end.insert(tup.0, idx);
        }
        if tup.2 < 0 || tup.2 == 1 {
            right_end.insert(tup.0, idx);
        }
    }

    let mut ret = Vec::new();

    for node in &graph.bnodes {
        ret.push((
            node.id,
            *left_end
                .get(&node.id)
                .expect("node is not present in left end list"),
            *right_end
                .get(&node.id)
                .expect("node is not present in right end list"),
        ));
    }

    return ret;
}

/// The iterative version of the main algorithm.
/// Let us recall the format of the input:
///     - ints: contains the nice interval representation, i.e.
///     elements of the form (y, a_y, b_y).
///     - crossing_dicts: given a tuple (u,v), tells you the number
///     of crossings between edges adjacent to u and edges adjacent to v
///     when placing u before v.
///     - k: the upper allowed limit to the number of crossings
///
/// The result is an Ok value if there is indeed an ordering allowing
/// less than k crossings, and Err otherwise.
pub fn kobayashi_tamaki(
    graph: &Graph,
    crossing_dict: &HashMap<(usize, usize), usize>,
) -> Result<Vec<usize>, String> {
    if graph.bnodes.len() == 1 {
        let mut vec = Vec::new();
        vec.push(graph.bnodes[0].id);
        return Ok(vec);
    }



    let ints = nice_interval_repr(graph);

    // 2|Y| in the article. we start at 0 here.
    let max_t: usize = 2 * ints.len();

    // computing the sets L_t and M_t, recording when t=a_y or t=b_y
    let mut m: Vec<Vec<usize>> = vec![Vec::new(); max_t]; // t to sorted list of Mt elements
    let mut l: HashMap<usize, Vec<usize>> = HashMap::new(); // t to sorted list of Lt elements
    for t in 0..max_t {
        l.entry(t).or_default();
    }
    let mut a: HashMap<usize, usize> = HashMap::new(); // t to y such that t=a_y
    let mut b: HashMap<usize, usize> = HashMap::new(); // t to y such that t=b_y
    for tup in &ints {
        a.insert(tup.1, tup.0);
        b.insert(tup.2, tup.0);
        for t in (tup.1)..(tup.2) {
            // or_default puts an empty vector if t entry does not exist.
            match m[t].binary_search(&tup.0) {
                Ok(_) => {}
                Err(position) => m[t].insert(position, tup.0),
            }
        }
        for t in (tup.2)..max_t {
            l.entry(t).or_default().push(tup.0)
        }
    }

    // Looking at size of largest M_t
    let mut mt_sizes = Vec::new();
    //    let mut h: u32 = 0;
    for t in 0..max_t {
        let v = &m[t]; //.get(&t).unwrap();
        mt_sizes.push(v.len());
        //        h = std::cmp::max(h, v.len() as u32);
    }

    // pos
    let mut pos: HashMap<(usize, usize), usize> = HashMap::new();
    for t in 0..max_t {
        for (p, u) in m[t].iter().enumerate() {
            pos.insert((t, *u), p);
        }
    }

    // lt crossings
    let mut lt_crossings: Vec<Vec<usize>> = vec![Vec::new(); max_t];
    for t in 0..max_t {
        for u in m[t].iter() {
            lt_crossings[t].push(set_crossing(&l[&t], *u, &crossing_dict))
        }
    }

    // pointer towards best assignment
    let mut ptr: HashMap<(usize, usize), usize> = HashMap::new();

    // the vec. playing the role of opt[t]
    let mut vec_b: Vec<usize> = vec![0; 1 << mt_sizes[0]];

    // first iteration
    let y = a
        .get(&0)
        .expect("the first thing to happen should be an interval opening");
    ptr.insert((0, 1), *y);

    // filling table
    for t in 1..(m.len()) {
        let vec_a = vec_b.clone();
        vec_b.resize(1 << mt_sizes[t], 0);

        // storing local scoring values to avoid calling HashMap 2**mt_sizes[t] times
        let mut index_wise_crossing: Vec<Vec<usize>> = vec![vec![0; m[t].len()]; m[t].len()];
        for p in 0..m[t].len() {
            for p2 in 0..m[t].len() {
                if p != p2 {
                    index_wise_crossing[p][p2] = crossing_dict[&(m[t][p], m[t][p2])];
                }
            }
        }

        // if t is b_y for some y
        if let Some(y) = b.get(&t) {
            let p = pos.get(&(t - 1, *y)).unwrap();

            for x in 0..(1 << mt_sizes[t]) {
                // x = ***#### (pos = 4)
                // tmp = ***
                // tmp2 = ####
                // new_x = ***1 then ***1#### (y has been inserted at its pos. it is counted
                // starting from the en)
                let tmp = x >> p; // ***
                let mask = tmp << p; // ***0000
                let tmp2 = x ^ mask; // 000####
                let mut new_x = (tmp << 1) + 1; // ***1
                new_x = (new_x << p) + tmp2; // ***10000+####

                vec_b[x] = vec_a[new_x];
            }
        }

        // if t is a_y for some y
        if let Some(y) = a.get(&t) {
            let p = pos.get(&(t, *y)).unwrap();

            for x in 0..(1 << mt_sizes[t]) {
                if (x >> p & 1) == 0 {
                    // if x does not contain y
                    // ###0**** --> ###**** (here p=4)
                    let mut tmp = x >> p + 1;
                    let mask = tmp << p + 1;
                    let tmp2 = x ^ mask;
                    tmp = tmp << p;
                    let new_x = tmp + tmp2;

                    vec_b[x] = vec_a[new_x];
                } else {
                    let mut best_sc = usize::MAX;
                    let mut best_x = usize::MAX;
                    for p in 0..mt_sizes[t] {
                        // if m[p] not in the set represented by x
                        if (x >> p) & 1 == 0 {
                            continue;
                        }
                        // ****1### -> ****0### (here p=3)
                        let mut tmp = x >> p; // ****1
                        let mask = tmp << p; // ****1000
                        let tmp2 = x ^ mask; // 00000###
                        tmp = (tmp >> 1) << 1; //****0 removing last bit.
                        let index = (tmp << p) + tmp2; //****0###

                        let mut sc: usize = vec_b[index];
                        sc += lt_crossings[t][p];
                        for p2 in 0..mt_sizes[t] {
                            if (x >> p2) & 1 == 0 {
                                continue;
                            }
                            if p != p2 {
                                //                                sc += crossing_dict[&(m[t][p2],m[t][p])]
                                sc += index_wise_crossing[p2][p];
                            }
                        }

                        if sc < best_sc {
                            best_sc = sc;
                            best_x = m[t][p];
                        }
                    }
                    vec_b[x] = best_sc;
                    ptr.insert((t, x), best_x);
                }
            }
        }
    }

    // reconstructing solution
    let mut solution = Vec::new();
    let mut s: Vec<usize> = Vec::new();
    let mut t: usize = max_t - 1;

    loop {
        if t == 0 && s.len() == 0 {
            break;
        }
        // if t is b_y for some y
        if let Some(y) = b.get(&t) {
            t -= 1;
            match s.binary_search(&y) {
                Ok(_pos) => {}
                Err(pos) => s.insert(pos, *y),
            }
        }

        // if t is a_y for some y
        if let Some(y) = a.get(&t) {
            if !s.contains(&y) {
                t -= 1;
            } else {
                if t == 0 {
                    assert_eq!(s.len(), 1);
                }
                let best_x = ptr.get(&(t, vec2int(&s, &m[t]))).unwrap();
                solution.push(*best_x);
                s.remove(s.binary_search(best_x).unwrap());
            }
        }
    }
    solution.reverse();

    Ok(solution) // It is a return (see "expressions" in rust)
}

/// The struct containing all the necessary "constant info"
/// that needs to be passed along the recursive calls
/// in a memoized version of the main algorithm
/// (Kobayashi-Tamaki). A reference will be passed,
/// and the recursive function will only need to read
/// from this object.
#[derive(Debug)]
struct ConstantInfo {
    b: HashMap<usize, usize>,
    a: HashMap<usize, usize>,
//    inv_a: HashMap<usize, usize>,
    pos: HashMap<(usize, usize), usize>,
    mt_sizes: Vec<usize>,
    lt_crossings: Vec<Vec<usize>>,
    crossing_dict: HashMap<(usize, usize), usize>,
    m: Vec<Vec<usize>>,
//    l: HashMap<usize, Vec<usize>>,
//    first_b: usize,
}

/// The function recomputing strongly connected components
/// in the sub-instance (t,x). It returns a integer representing
/// a subset of x, which are "candidates". When looking for the
/// last element in an optimal ordering of the sub-instance
/// associated to (t,x), it is enough to only look
/// at candidates.
//fn recomputing_scc_rule(t: usize, x: usize, constant_info: &ConstantInfo) -> usize {
//    let mut h = petgraph::graph::Graph::<usize, usize>::new();
//
//    // mapping node in h to u balue
//    let mut map: HashMap<usize, usize> = HashMap::new();
//    let mut inv_map: HashMap<usize, petgraph::graph::NodeIndex> = HashMap::new();
//
//    if constant_info.l[&t].len() > 0 {
//        // node representing all lt
//        let nu = h.add_node(usize::MAX);
//        map.insert(nu.index(), usize::MAX);
//        inv_map.insert(usize::MAX, nu);
//    }
//
//    // loop over nodes of m[t] to add edges with usize::UMAX (lt nodes)
//    for p in 0..constant_info.mt_sizes[t] {
//        // if m[p] not in the set represented by x
//        if (x >> p) & 1 == 0 {
//            continue;
//        }
//        let u = constant_info.m[t][p];
//        let nu = h.add_node(u);
//        map.insert(nu.index(), u);
//        inv_map.insert(u, nu);
//
//        for ul in &constant_info.l[&t] {
//            if let Some(c) = constant_info.crossing_dict.get(&(u, *ul)) {
//                let c2 = constant_info.crossing_dict.get(&(*ul, u)).unwrap();
//                if c < c2 {
//                    let a = *inv_map.get(&usize::MAX).unwrap();
//                    let b = *inv_map.get(&u).unwrap();
//                    if !h.contains_edge(b, a) {
//                        h.add_edge(b, a, 1);
//                        break;
//                    }
//                }
//                if c > c2 {
//                    let a = *inv_map.get(&usize::MAX).unwrap();
//                    let b = *inv_map.get(&u).unwrap();
//                    if !h.contains_edge(a, b) {
//                        h.add_edge(a, b, 1);
//                        break;
//                    }
//                }
//            }
//        }
//        if *constant_info.inv_a.get(&u).unwrap() > constant_info.first_b {
//            let a = *inv_map.get(&usize::MAX).unwrap();
//            let b = *inv_map.get(&u).unwrap();
//            if !h.contains_edge(a, b) {
//                h.add_edge(a, b, 1);
//            }
//        }
//    }
//
//    // adding edges inside m[t]
//    for p in 0..constant_info.mt_sizes[t] {
//        if (x >> p) & 1 == 0 {
//            continue;
//        }
//        for p2 in 0..constant_info.mt_sizes[t] {
//            if (x >> p2) & 1 == 0 {
//                continue;
//            }
//            if p != p2 {
//                let u = constant_info.m[t][p];
//                let v = constant_info.m[t][p2];
//                let a = *inv_map.get(&u).unwrap();
//                let b = *inv_map.get(&v).unwrap();
//                let c = constant_info.crossing_dict.get(&(u, v));
//                let c2 = constant_info.crossing_dict.get(&(v, u));
//                if c < c2 {
//                    if !h.contains_edge(a, b) {
//                        h.add_edge(a, b, 1);
//                    }
//                }
//                if c > c2 {
//                    if !h.contains_edge(b, a) {
//                        h.add_edge(b, a, 1);
//                    }
//                }
//            }
//        }
//    }
//
//    let sccs = tarjan_scc(&h);
//    // then the last scc in topo order (first in sccs) is last
//
//    let mut candidates = 0;
//
//    // putting vertices in graph
//    for u in sccs[0].clone().into_iter() {
//        let id = map.get(&u.index()).unwrap();
//        if *id != usize::MAX {
//            let p = constant_info.pos.get(&(t, *id)).unwrap();
//            candidates += 1 << p;
//        }
//    }
//
//    return candidates;
//}

/// bit manipulation setting the p-th bit of x to 0.
fn turn_pth_bit_off(x: usize, p: usize) -> usize {
    // ****1### -> ****0### (here p=3)
    let mut tmp = x >> p; // ****1
    let mask = tmp << p; // ****1000
    let tmp2 = x ^ mask; // 00000###
    tmp = (tmp >> 1) << 1; //****0 removing last bit.
    let index = (tmp << p) + tmp2; //****0###
    return index;
}

/// Memoized recursive function for the computation
/// of the optimal number of crossings using Kobayashi-Tamaki.
fn opt_num_crossings(
    t: usize,
    x: usize,
    dp_table: &mut HashMap<(usize, usize), usize>,
    ptr: &mut HashMap<(usize, usize), Vec<usize>>,
    constant_info: &ConstantInfo,
) -> usize {
    //    info!("calling on {:?}", (t,x));
    if let Some(value) = dp_table.get(&(t, x)) {
        return *value;
    }

    if t == 0 && x == 0 {
        dp_table.insert((t, x), 0);
        return *dp_table.get(&(t, x)).unwrap();
    }

    // if t is b_y for some y
    if let Some(y) = constant_info.b.get(&t) {
        let p = constant_info.pos.get(&(t - 1, *y)).unwrap();

        let tmp = x >> p; // ***
        let mask = tmp << p; // ***0000
        let tmp2 = x ^ mask; // 000####
        let mut new_x = (tmp << 1) + 1; // ***1
        new_x = (new_x << p) + tmp2; // ***10000+####

        let opt = opt_num_crossings(t - 1, new_x, dp_table, ptr, constant_info);
        dp_table.insert((t, x), opt);
    }

    // if t is a_y for some y
    if let Some(y) = constant_info.a.get(&t) {
        let p = constant_info.pos.get(&(t, *y)).unwrap();

        if (x >> p & 1) == 0 {
            // if x does not contain y
            // ###0**** --> ###**** (here p=4)
            let mut tmp = x >> p + 1;
            let mask = tmp << p + 1;
            let tmp2 = x ^ mask;
            tmp = tmp << p;
            let new_x = tmp + tmp2;
            let opt = opt_num_crossings(t - 1, new_x, dp_table, ptr, constant_info);
            dp_table.insert((t, x), opt);
        } else {
//            // RECOMPUTING SCC RULE
//            let candidates = recomputing_scc_rule(t, x, constant_info);
//            // END RECOMPUTING SCC RULE
            let candidates = x;

            let mut best_sc = usize::MAX;
            let mut best_x = usize::MAX;
            for p in 0..constant_info.mt_sizes[t] {
                // if m[p] not in the set represented by x
                if (x >> p) & 1 == 0 {
                    continue;
                }
                // if m[p] is not a candidate
                if (candidates >> p) & 1 == 0 {
                    continue;
                }
                // ****1### -> ****0### (here p=3)
                let index = turn_pth_bit_off(x, p);

                let mut sc: usize = opt_num_crossings(t, index, dp_table, ptr, constant_info);
                sc += constant_info.lt_crossings[t][p];
                for p2 in 0..constant_info.mt_sizes[t] {
                    if (x >> p2) & 1 == 0 {
                        continue;
                    }
                    if p != p2 {
                        sc += constant_info.crossing_dict
                            [&(constant_info.m[t][p2], constant_info.m[t][p])];
                    }
                }

                if sc < best_sc {
                    best_sc = sc;
                    best_x = constant_info.m[t][p];
                }
            }
            dp_table.insert((t, x), best_sc);
            ptr.insert((t, x), vec![best_x; 1]);
        }
    }

    return *dp_table.get(&(t, x)).unwrap();
}

/// Backtrace of the memoization implementation of
/// Kobayashi-Tamaki, which reconstructs the
/// final solution
fn backtrace(
    t: usize,
    x: usize,
    ptr: &HashMap<(usize, usize), Vec<usize>>,
    constant_info: &ConstantInfo,
) -> Vec<usize> {
    if t == 0 && x == 0 {
        return Vec::new();
    }

    let mut return_value: Vec<usize> = Vec::new();

    // if t is b_y for some y
    if let Some(y) = constant_info.b.get(&t) {
        let p = constant_info.pos.get(&(t - 1, *y)).unwrap();

        let tmp = x >> p; // ***
        let mask = tmp << p; // ***0000
        let tmp2 = x ^ mask; // 000####
        let mut new_x = (tmp << 1) + 1; // ***1
        new_x = (new_x << p) + tmp2; // ***10000+####

        // return to lighten up table a bit
        return_value = backtrace(t - 1, new_x, ptr, constant_info);
    }

    // if t is a_y for some y
    if let Some(y) = constant_info.a.get(&t) {
        let p = constant_info.pos.get(&(t, *y)).unwrap();

        if (x >> p & 1) == 0 {
            // if x does not contain y
            // ###0**** --> ###**** (here p=4)
            let mut tmp = x >> p + 1;
            let mask = tmp << p + 1;
            let tmp2 = x ^ mask;
            tmp = tmp << p;
            let new_x = tmp + tmp2;
            return_value = backtrace(t - 1, new_x, ptr, constant_info);
        } else {
            let extension: Vec<usize> = ptr.get(&(t, x)).unwrap().to_vec();
            let mut new_set: Vec<usize> = Vec::new();
            for elem in int2vec(x, &constant_info.m[t]) {
                if !extension.contains(&elem) {
                    new_set.push(elem);
                }
            }
            let new_index: usize = vec2int(&new_set, &constant_info.m[t]);
            let mut rest: Vec<usize> = backtrace(t, new_index, ptr, constant_info);
            rest.extend_from_slice(extension.as_slice());
            return_value = rest;
        }
    }
    return return_value;
}

pub fn total_instance_size(graph: &Graph) -> usize {
    
    let ints = nice_interval_repr(graph);

    // 2|Y| in the article. we start at 0 here.
    let max_t: usize = 2 * ints.len();

    // computing the sets L_t and M_t, recording when t=a_y or t=b_y
    let mut m: Vec<Vec<usize>> = vec![Vec::new(); max_t]; // t to sorted list of Mt elements
                                                          //
    for tup in &ints {
        for t in (tup.1)..(tup.2) {
            // or_default puts an empty vector if t entry does not exist.
            match m[t].binary_search(&tup.0) {
                Ok(_) => {}
                Err(position) => m[t].insert(position, tup.0),
            }
        }
    }
    
    // Looking at size of largest M_t
    let mut mt_sizes = Vec::new();
    //    let mut h: u32 = 0;
    for t in 0..max_t {
        let v = &m[t]; //.get(&t).unwrap();
        mt_sizes.push(v.len());
        //        h = std::cmp::max(h, v.len() as u32);
    }

    let mut total_instance_size = 0;
    for s in &mt_sizes {
        total_instance_size += 1 << s;
    }
    return total_instance_size;

}

/// memoization version of Kobayashi-Tamaki.
/// The preambule is the same as the iterative
/// version.
pub fn recursive_kt(
    graph: &Graph,
    crossing_dict: &HashMap<(usize, usize), usize>,
) -> Result<Vec<usize>, String> {
    if graph.bnodes.len() == 1 {
        let mut vec = Vec::new();
        vec.push(graph.bnodes[0].id);
        return Ok(vec);
    }

    let ints = nice_interval_repr(graph);

    // 2|Y| in the article. we start at 0 here.
    let max_t: usize = 2 * ints.len();

    // computing the sets L_t and M_t, recording when t=a_y or t=b_y
    let mut m: Vec<Vec<usize>> = vec![Vec::new(); max_t]; // t to sorted list of Mt elements
    let mut l: HashMap<usize, Vec<usize>> = HashMap::new(); // t to sorted list of Lt elements
    for t in 0..max_t {
        l.entry(t).or_default();
    }
    let mut a: HashMap<usize, usize> = HashMap::new(); // t to y such that t=a_y
    let mut inv_a: HashMap<usize, usize> = HashMap::new(); // y to t such that t=a_y
    let mut b: HashMap<usize, usize> = HashMap::new(); // t to y such that t=b_y
    for tup in &ints {
        a.insert(tup.1, tup.0);
        inv_a.insert(tup.0, tup.1);
        b.insert(tup.2, tup.0);
        for t in (tup.1)..(tup.2) {
            // or_default puts an empty vector if t entry does not exist.
            match m[t].binary_search(&tup.0) {
                Ok(_) => {}
                Err(position) => m[t].insert(position, tup.0),
            }
        }
        for t in (tup.2)..max_t {
            l.entry(t).or_default().push(tup.0)
        }
    }

    // Looking at size of largest M_t
    let mut mt_sizes = Vec::new();
    //    let mut h: u32 = 0;
    for t in 0..max_t {
        let v = &m[t]; //.get(&t).unwrap();
        mt_sizes.push(v.len());
        //        h = std::cmp::max(h, v.len() as u32);
    }

    // pos
    let mut pos: HashMap<(usize, usize), usize> = HashMap::new();
    for t in 0..max_t {
        for (p, u) in m[t].iter().enumerate() {
            pos.insert((t, *u), p);
        }
    }

    // lt crossings
    let mut lt_crossings: Vec<Vec<usize>> = vec![Vec::new(); max_t];
    for t in 0..max_t {
        for u in m[t].iter() {
            lt_crossings[t].push(set_crossing(&l[&t], *u, &crossing_dict))
        }
    }

    // first b
    let mut first_b: usize = usize::MAX;
    for b in b.keys() {
        first_b = std::cmp::min(first_b, *b);
    }

    // END PREAMBULE
    //

    let constant_info = ConstantInfo {
        b: b,
        a: a,
//        inv_a: inv_a,
        pos: pos,
        mt_sizes: mt_sizes,
        lt_crossings: lt_crossings,
        crossing_dict: crossing_dict.clone(),
        m: m,
//        l: l,
//        first_b: first_b,
    };

    // init dp table
    let mut dp_table: HashMap<(usize, usize), usize> = HashMap::new();
    let mut ptr: HashMap<(usize, usize), Vec<usize>> = HashMap::new();

    // computing optimal number of crossings
    let _opt = opt_num_crossings(max_t - 1, 0, &mut dp_table, &mut ptr, &constant_info);

    // reconstructing optimal solution
    let ordering = backtrace(max_t - 1, 0, &ptr, &constant_info);

    return Ok(ordering);
}



































//outs MUST contain an entry for each i in 0..nbnodes
pub fn write_digraph( nbnodes: usize, outs: &HashMap<usize, HashSet<usize>>, 
                      arc_weights : &HashMap< (usize, usize), usize>, filename : &String )
{
    let mut data_file = File::create(filename).expect("Failed to create file");
    
    //line 1 = nb vertices
    let _ = write!(data_file, "{}\n", nbnodes);
    
    //other lines = arcs, format "v1 v2 weight", where v1 v2 are vertex numbers, indexed at 1
    for i in 0 .. nbnodes{
        for j in &outs[&i]{
            let _ = write!(data_file, "{} {} {}\n", i + 1, j + 1, arc_weights[ &(i, *j) ]);
        }
    }
}








//taken from the web
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}


pub fn solve_dfas_from_file( filename : &String ) -> ( usize, Vec<usize> )
{
    let mut infofilename : String = String::new();
    infofilename += filename;
    infofilename += ".info";
    
    info!("Calling ./dfas_v2/dfas --in={} --out={}", filename, infofilename);
    
    let _output = Command::new("./dfas_v2/dfas")
                     .arg("--in=".to_owned() + &filename)
                     .arg("--out=".to_owned() + &infofilename)
                     .output()
                     .expect("failed to execute process");
                    

    let mut cost : usize = 0;
    let mut topo_sort : Vec<usize> = Vec::new();
    
    if let Ok(lines) = read_lines(infofilename) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines.flatten() {
            
            //info!("Line = {}", line);
            
            let parts : Vec<&str> = line.split("=").collect();
            
            if parts.len() >= 2
            {
                if parts[0] == "cost"
                {
                    
                    cost = parts[1].parse::<usize>().unwrap();
                }
                else if parts[0] == "topo_sort"
                {
                    let topo_vertices : Vec<&str> = parts[1].split(" ").collect();
                    
                    for t in topo_vertices
                    {
                        if t != ""
                        {
                            let v : usize = t.parse::<usize>().unwrap();
                            topo_sort.push(v);
                        }
                    }
                }
            }
            
        }
    }
    
    return (cost, topo_sort);
}




















pub fn solve_dfas( graph: &Graph, crossing_dict: &HashMap<(usize, usize), usize>, instance_id : usize ) -> (usize, Vec<usize>)
{
    //we are essentially rebuilding the digraph from compute_scc, but whatever...
    //because maxsatsolver wants var names 1,2,...,n, we work with array indices of graph.bnodes instead of id values
    
    
    //compute map of out-neighbors (adjacency list with maps).  
    //key = index in array of graph.bnodes, value = set of array indices of out-neighbors
    let mut outs: HashMap<usize, HashSet<usize>> = HashMap::new();
    let mut arc_weights: HashMap< (usize, usize), usize > = HashMap::new();
    
    
    for i in 0 .. graph.bnodes.len()
    {
        let i_id = graph.bnodes[i].id;
        
        let mut i_outs : HashSet<usize> = HashSet::new();
        
        
        for j in 0 .. graph.bnodes.len()
        {
            let j_id = graph.bnodes[j].id;
            
            if let Some(c_ji) = crossing_dict.get(&(j_id, i_id)) {
                let c_ij = crossing_dict.get(&(i_id, j_id)).unwrap();
                if *c_ji > *c_ij {
                    //j < i has more crossings than i < j => arc (i, j)
                    i_outs.insert(j);
                    arc_weights.insert( (i, j), *c_ji - *c_ij );
                }
            }
            
            if !i_outs.contains(&j){
                //check same condition as compute_scc for when i is forced left of j
                let u = &graph.bnodes[i];
                let v = &graph.bnodes[j];
                if u.right <= v.left && u.left < v.right {
                    i_outs.insert(j);
                    arc_weights.insert( (i, j), 30000);
                }
            }
            
        }
        
        outs.insert(i, i_outs);
        
    }
    
    
    //create tmp dir if not exists
    let _ = fs::create_dir_all("./tmp");
    
    let mut filename : String = String::new();
    filename += "tmp/scc";
    filename +=  instance_id.to_string().as_str();
    filename += ".gr";
    
    write_digraph( graph.bnodes.len(), &outs, &arc_weights, &filename );
    
    let (cost, topo_sort) = solve_dfas_from_file( &filename );
    
    
    let mut toposort_bnode_ids : Vec<usize> = Vec::new();
    
    for t in &topo_sort{
        toposort_bnode_ids.push( graph.bnodes[*t].id );
    }
    
    
    //This checks that the claimed cost is the same as crossing_dict would calculate
    let mut cost_lb : usize = 0;
    let mut cost_maxsat : usize = 0;
    for i in 0 .. toposort_bnode_ids.len()
    {
        let id_i = toposort_bnode_ids[i];
        let bnode_i = &graph.bnodes[ topo_sort[i] ];
        
        //for each j before i in topo sort
        for j in 0 .. i 
        {
            let id_j = toposort_bnode_ids[j];
            let bnode_j = &graph.bnodes[ topo_sort[j] ];
            
            //case 1 : i is completely left of j, but got placed after -> should never happen since arc weights are super high
            //also contribution to lb is 0
            //condition is copied from h construction
            if bnode_i.right <= bnode_j.left && bnode_i.left < bnode_j.right
            {
                info!("Error: i.right is <= j.left but i after j, i={}  j={}", i, j);
            }
            
            
            if crossing_dict.contains_key(&(id_i, id_j))
            {
                let cost_ij = crossing_dict.get(&(id_i, id_j)).unwrap();
                let cost_ji = crossing_dict.get(&(id_j, id_i)).unwrap();
                
                cost_lb += std::cmp::min(cost_ij, cost_ji);            
                cost_maxsat += cost_ji;    
            }
            
        }
    }
    
    info!("total cost lb={}  maxsat={}  diff={}", cost_lb, cost_maxsat, cost_maxsat - cost_lb);
    
    if cost_maxsat - cost_lb != cost
    {
        info!("ERROR: cost_maxsat - cost_lb != cost");
    }
    
    
    
    return (cost, toposort_bnode_ids)
    
}





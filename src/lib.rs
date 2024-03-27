use clap::Parser;
use petgraph::algo::*;
use std::collections::HashMap;
use std::collections::HashSet;
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

pub fn find_twins(graph: &Graph) -> HashMap<usize, usize> {
    // main hash table
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

pub fn compute_scc(graph: &Graph, crossing_dict: &HashMap<(usize, usize), usize>) -> Vec<Graph> {
    // just mapping id to bnode to remember it
    let mut id_to_bnodes: HashMap<usize, &BNode> = HashMap::new();
    for u in &graph.bnodes {
        id_to_bnodes.insert(u.id, u);
    }
    //    println!("id_to_bnodes {:?}", id_to_bnodes);

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
                    h.add_edge(a, b, *c-*c2);
                }
            } else if *c < *c2 {
                let a = *map.get(v).unwrap();
                let b = *map.get(u).unwrap();
                if !h.contains_edge(a, b) {
                    h.add_edge(a, b, *c2-*c);
                }
            }
        }
    }
    //    println!("simplified graph");
    //    println!("{:?}", petgraph::dot::Dot::with_config(&h, &[petgraph::dot::Config::NodeIndexLabel]));
    //    println!("end simplified graph");

    for u in &graph.bnodes {
        for v in &graph.bnodes {
            if u.id != v.id {
                if u.right <= v.left && u.left < v.right {
                    let a = *map.get(&u.id).unwrap();
                    let b = *map.get(&v.id).unwrap();
                    if !h.contains_edge(a, b) {
                        h.add_edge(a, b, 0);
                    }
                }
            }
        }
    }

    //    for u in &graph.bnodes {
    //        for v in &graph.bnodes {
    //            if u.id !=v.id {
    //                let a = *map.get(&u.id).unwrap();
    //                let b = *map.get(&v.id).unwrap();
    //                if !h.contains_edge(a,b) && !h.contains_edge(b,a) {
    //                    assert_eq!(*crossing_dict.get(&(u.id, v.id)).unwrap(),*crossing_dict.get(&(v.id, u.id)).unwrap());
    //                }
    //            }
    //        }
    //    }
    println!(
        "h has {:?} edges over {:?} possible",
        h.edge_count(),
        h.node_count() * (h.node_count() - 1) / 2
    );

    //    println!("map {:?}", map);
    //    println!("directed graph {:?}", h);
    let sccs = tarjan_scc(&h);
    //    println!("sccs {:?}", sccs);
    //println!("{:?}", petgraph::dot::Dot::with_config(&h, &[petgraph::dot::Config::GraphContentOnly]));

    let mut graph_vec: Vec<Graph> = Vec::new();
    for scc in &sccs {
        // dot visu
                let mut scc_h = h.clone();
                scc_h.retain_nodes(|scc_h,u| scc.contains(&u));
         println!("{:?}", petgraph::dot::Dot::with_config(&scc_h, &[petgraph::dot::Config::GraphContentOnly]));

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

//fn int2vec(x: usize, v: &Vec<usize>) -> Vec<usize> {
//    let mut s = Vec::new();
//    for (i, y) in v.iter().enumerate() {
//        // if i-th bit of x is 1, y in s
//        if (x >> i) & 1 != 0 {
//            s.push(*y);
//        }
//    }
//    return s;
//}

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

/// The main algorithm. Let us recall the format of the input:
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

    let mut total_number_operations = 0;
    for t in 1..(m.len()) {
        total_number_operations +=  mt_sizes[t];
    }

    let mut ops = 0;
    let mut thresh = 0.01;
    let inc = 0.01;

    println!("mt sizes {:?}", mt_sizes);

    // filling table
    for t in 1..(m.len()) {

        // estimating fraction of computation left
        ops += mt_sizes[t-1];
        println!("{:?} out of {:?}", ops, total_number_operations);
//        if ops as f64 > thresh*(total_number_operations as f64) {
//            print!("...{:?}", (thresh*100.0).floor());
//            thresh += inc;
//        }

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
    println!("");

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

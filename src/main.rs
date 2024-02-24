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
                println!("solution {:?}", vec);
                // Writing result in output file (name same as input, extension changed)
                let v: Vec<String> = vec.into_iter().map(|x| x.to_string()).collect();
                let _ = std::fs::write(outname, v.join("\n"));
                break;
            }
            Err(_) => {
                // if not possible with continue to next loop to try k+1
                k += 1;
                continue;
            }
        }
    }
}

/// Computing the number of binary numbers having Hamming
/// weight w and L bits in their decomposition (with potential
/// 0s on the left side)
//fn num_fixed_hamming_weigth(L: usize, w: usize) -> u64 {
//
//}
//
//fn subset_rank(subset: Vec<bool>) -> u64 {
//
//
//}
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

fn vec2int(s: &Vec<usize>, v: &Vec<usize>) -> usize {
    let mut m: HashMap<usize, usize> = HashMap::new();
    for (i, y) in v.iter().enumerate() {
        m.insert(*y, i);
    }
    let mut x: usize = 0;
    for y in s {
        x += usize::pow(2, m.remove(&y).unwrap().try_into().unwrap());
    }
    return x;
}

fn set_crossing(s: &Vec<usize>, x: usize, crossing_dict: &HashMap<(usize, usize), u32>) -> u32 {
    let mut c: u32 = 0;
    for u in s {
        c += crossing_dict[&(*u, x)];
    }
    return c;
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
fn kobayashi_tamaki(
    ints: &Vec<(usize, usize, usize)>,
    crossing_dict: &HashMap<(usize, usize), u32>,
    k: u32,
) -> Result<Vec<usize>, String> {
    // 2|Y| in the article. we start at 0 here.
    let max_t: usize = 2 * ints.len();

    // computing the sets L_t and M_t, recording when t=a_y or t=b_y
    let mut m: HashMap<usize, Vec<usize>> = HashMap::new(); // t to sorted list of Mt elements
    let mut l: HashMap<usize, Vec<usize>> = HashMap::new(); // t to sorted list of Lt elements
    for t in 0..max_t {
        l.entry(t).or_default();
        m.entry(t).or_default();
    }
    let mut a: HashMap<usize, usize> = HashMap::new(); // t to y such that t=a_y
    let mut b: HashMap<usize, usize> = HashMap::new(); // t to y such that t=b_y
    for tup in ints {
        a.insert(tup.1, tup.0);
        b.insert(tup.2, tup.0);
        for t in (tup.1)..(tup.2) {
            // or_default puts an empty vector if t entry does not exist.
            m.entry(t).or_default().push(tup.0);
        }
        for t in (tup.2)..max_t {
            l.entry(t).or_default().push(tup.0)
        }
    }

    println!("M: {:?}",m);
    // Looking at size of largest M_t
    let mut mt_sizes = Vec::new();
    let mut h: u32 = 0;
    for t in 0..max_t {
        let v = m.get(&t).unwrap();
        mt_sizes.push(v.len() as u32);
        h = std::cmp::max(h, v.len() as u32);
    }

    // This is where k comes into account
    if h * (h - 1) > k {
        return Err("k not high enough".to_string());
    }

    println!("M sets {:?}", m);
    println!("L sets {:?}", l);

    // dynamic programming table initialization
    let mut opt = vec![Vec::new(); m.len()];
    for (t, v) in opt.iter_mut().enumerate() {
        let new_size = usize::pow(2, mt_sizes[t]); // size depends on t
        v.resize(new_size, 0); // init at zero (!)
    }


    // filling table
    for t in 1..(m.len()) {
        // TODO: iteration below is naive, and inefficient in terms of contiguity.
        for x in 0..opt[t].len() {
            let mut s = int2vec(x, &m[&t]);

            // if t is b_y for some y
            if let Some(y) = b.get(&t) {
                match s.binary_search(&y) {
                    Ok(_) => {}
                    Err(pos) => s.insert(pos, *y),
                }
                opt[t][x] = opt[t - 1][vec2int(&s, &m[&(t - 1)])]
            }

            // if t is a_y for some y
            if let Some(y) = a.get(&t) {
                if !s.contains(&y) {
                    opt[t][x] = opt[t - 1][vec2int(&int2vec(x, &m[&t]), &m[&(t - 1)])];
                } else {
                    let mut best_sc = u32::MAX;
                    for x in &s {
                        let mut s2 = s.clone();
                        s2.remove(s2.binary_search(&x).unwrap());

                        let mut sc: u32 = opt[t][vec2int(&s2, &m[&t])];
                        sc += set_crossing(&l[&t], *x, &crossing_dict);
                        sc += set_crossing(&s2, *x, &crossing_dict);

                        if sc < best_sc {
                            best_sc = sc;
                        }
                    }
                    opt[t][x] = best_sc;
                }
            }
        }
    }
    println!("table {:?}", opt);
    println!("OPT {:?}", opt[max_t-1][0]);
    if opt[max_t-1][0] > k {
        return Err("k not high enough".to_string());
    }

    // reconstructing solution
    let mut solution = Vec::new();
    let mut s: Vec<usize> = Vec::new();
    let mut t: usize = max_t-1;

    loop {
        if t==0 && s.len()==0 {
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
                let mut best_sc = u32::MAX;
                let mut best_x = usize::MAX;
                for x in &s {
                    let mut s2 = s.clone();
                    s2.remove(s2.binary_search(&x).unwrap());

                    let mut sc: u32 = opt[t][vec2int(&s2, &m[&t])];
                    sc += set_crossing(&l[&t], *x, &crossing_dict);
                    sc += set_crossing(&s2, *x, &crossing_dict);

                    if sc < best_sc {
                        best_sc = sc;
                        best_x = *x;
                    }
                }
                solution.push(best_x);
                s.remove(s.binary_search(&best_x).unwrap());
            }
        }
    }
    solution.reverse();

    //    let mut v = Vec::new();
    //    for p in ints {
    //        v.push(p.0);
    //    }
    Ok(solution) // It is a return (see "expressions" in rust)
}

/// Computation of all crossing values c(u,v) for u,v two
/// vertices of Y, the side of the graph that must be ordered.
/// The result is an Hashmap associating keys (u,v) to
/// the integer value c(u,v), the number of crossings
/// obtained between edges adjacent to u and edges
/// adjacent to v if u is placed before v in the order.
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

/// a simple function that excludes the first component
/// of a tuple, so that it may be ignored in a lexico-graphic
/// ordering
fn exclude_first(p: &(usize, i32, i32)) -> (i32, i32) {
    (p.1, p.2)
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
        let pl: (usize, i32, i32) = (node.id, node.left as i32, 2*degree);
        let pr: (usize, i32, i32) = (node.id, node.right as i32, - degree);
        p.push(pl);
        p.push(pr);
    }

    // lexico-graphic order while ignoring first
    p.sort_by(|a, b| exclude_first(a).cmp(&exclude_first(b)));

    println!("sorted p: {:?}", p);

    let mut left_end = HashMap::new();
    let mut right_end = HashMap::new();

    for (idx, tup) in p.iter().enumerate() {
        if tup.2 > 0 {
            left_end.insert(tup.0, idx);
        }
        if tup.2 < 0 {
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

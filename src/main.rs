use clap::Parser;
use sherby_pace_2024::*;
use std::collections::HashMap;

fn main() {
    let args = Cli::parse();

    let graph: Graph = parse_graph(&args.graph);

    let ints = nice_interval_repr(&graph);

    // TODO: implement more efficient crossing values procedure
    let crossing_dict = crossing_values(&graph);

    println!("crossing values {:?}", crossing_dict);

    let mut outname = args.graph.clone();
    outname.set_extension("sol");

    // main call
    let vec = kobayashi_tamaki(&ints, &crossing_dict).unwrap();

    println!("solution {:?}", vec);
    // Writing result in output file (name same as input, extension changed)
    let v: Vec<String> = vec.into_iter().map(|x| x.to_string()).collect();
    let _ = std::fs::write(outname, v.join("\n"));
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

fn set_crossing(s: &Vec<usize>, x: usize, crossing_dict: &HashMap<(usize, usize), usize>) -> usize {
    let mut c: usize = 0;
    for u in s {
        c += crossing_dict.get(&(*u, x)).unwrap_or(&(0 as usize));
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
    crossing_dict: &HashMap<(usize, usize), usize>,
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

    println!("M: {:?}", m);
    // Looking at size of largest M_t
    let mut mt_sizes = Vec::new();
    let mut h: u32 = 0;
    for t in 0..max_t {
        let v = m.get(&t).unwrap();
        mt_sizes.push(v.len() as u32);
        h = std::cmp::max(h, v.len() as u32);
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
                    let mut best_sc = usize::MAX;
                    for x in &s {
                        let mut s2 = s.clone();
                        s2.remove(s2.binary_search(&x).unwrap());

                        let mut sc: usize = opt[t][vec2int(&s2, &m[&t])];
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
    println!("OPT {:?}", opt[max_t - 1][0]);

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
                let mut best_sc = usize::MAX;
                let mut best_x = usize::MAX;
                for x in &s {
                    let mut s2 = s.clone();
                    s2.remove(s2.binary_search(&x).unwrap());

                    let mut sc: usize = opt[t][vec2int(&s2, &m[&t])];
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
fn crossing_values(graph: &Graph) -> HashMap<(usize, usize), usize> {
    
    let mut d_less_than_x: HashMap<(usize,usize), usize> = HashMap::new();

    // #1 compute all d^<x(u) O(|X||Y|)
    for u in &graph.bnodes {
        let mut d: usize = 0;
        for x in &graph.anodes {
            d_less_than_x.insert((u.id,x.id),d);
            if u.neighbors.contains(&x.id) {
                d += 1;
            }
        }
    }

    // #2 compute all crossing values
    let mut crossing_dict: HashMap<(usize,usize), usize> = HashMap::new();
   
    for u in &graph.bnodes {
        for v in &graph.bnodes {
            if u.id==v.id {
                continue;
            }
            for x in &u.neighbors {
                if let Some(c) = crossing_dict.get(&(u.id,v.id)) {
                    crossing_dict.insert((u.id,v.id),*c+d_less_than_x.get(&(v.id,*x)).unwrap());
                }
                else {
                    crossing_dict.insert((u.id,v.id),*d_less_than_x.get(&(v.id,*x)).unwrap()); 
                }
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

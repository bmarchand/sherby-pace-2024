use clap::Parser;
use sherby_pace_2024::*;
use std::collections::HashMap;
use std::collections::HashSet;
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    let args = Cli::parse();

    let graph: Graph = parse_graph(&args.graph);
    
    let mut crossing_dict = orientable_crossing_values(&graph);

    let twin_mapping = find_twins(&graph);
    println!("{:?} twins out of {:?} vertices", twin_mapping.len(), graph.bnodes.len());

    let graph = process_twins_graph(graph, &twin_mapping);
    println!("num vertices after {:?}", graph.bnodes.len());

    crossing_dict = process_twins_crossing_dict(&twin_mapping, &crossing_dict);

    let ints = nice_interval_repr(&graph);

    let mut outname = args.graph.clone();
    outname.set_extension("sol");

    // main call
    let vec = kobayashi_tamaki(&ints, &crossing_dict).unwrap();

    let vec = add_twins(vec, &twin_mapping);

    // Writing result in output file (name same as input, extension changed)
    let v: Vec<String> = vec.into_iter().map(|x| x.to_string()).collect();
    let _ = std::fs::write(outname, v.join("\n"));

    let peak_mem = PEAK_ALLOC.peak_usage_as_mb();
    println!("peak memory: {} mb", peak_mem);
}


fn add_twins(vec: Vec<usize>, twin_mapping: &HashMap<usize,usize>) -> Vec<usize> {

    let mut inverse_mapping: HashMap<usize, Vec<usize>> = HashMap::new();

    for (u,v) in twin_mapping.iter() {
        if let Some(vector) = inverse_mapping.get_mut(v) {
            (*vector).push(*u);
        }
        else {
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


//fn process_twins_ints(twin_mapping: &HashMap<usize,usize>,ints: Vec<(usize, usize, usize)>) -> Vec<(usize,usize,usize)> {
//    
//    let mut new_ints: Vec<(usize,usize,usize)> = Vec::new();
//
//    for tup in &ints {
//        if twin_mapping.contains_key(&tup.0) {
//            continue
//        }
//        new_ints.push(*tup);
//    }
//    return new_ints;
//}

fn process_twins_crossing_dict(twin_mapping: &HashMap<usize,usize>, 
                 crossing_dict: &HashMap<(usize, usize), usize>) ->  HashMap<(usize, usize), usize>
{

    let mut new_crossing_dict: HashMap<(usize,usize), usize> = HashMap::new();

    for ((u,v), c) in crossing_dict.iter() {
        let mut w = u;
        let mut x = v;
        if let Some(y) = twin_mapping.get(u) {
            w = y;
        }
        if let Some(z) = twin_mapping.get(v) {
            x = z;
        }
        if let Some(c2) = new_crossing_dict.get(&(*w,*x)) {
            new_crossing_dict.insert((*w,*x), *c2+*c);
        }
        else {
            new_crossing_dict.insert((*w,*x), *c);
        }
    }
    return new_crossing_dict;
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
    if ints.len()==1 {
        let mut vec = Vec::new();
        vec.push(ints[0].0);
        return Ok(vec);
    }


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

    // Looking at size of largest M_t
    let mut mt_sizes = Vec::new();
    let mut h: u32 = 0;
    for t in 0..max_t {
        let v = m.get(&t).unwrap();
        mt_sizes.push(v.len() as u32);
        h = std::cmp::max(h, v.len() as u32);
    }

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

///// Computation of all crossing values c(u,v) for u,v two
///// vertices of Y, the side of the graph that must be ordered.
///// The result is an Hashmap associating keys (u,v) to
///// the integer value c(u,v), the number of crossings
///// obtained between edges adjacent to u and edges
///// adjacent to v if u is placed before v in the order.
/////
///// It works by first compuring $d^{<x}(u)$, the number of neighbors
///// of $u$ (a vertex of B) that are strictly before x, a vertex in A.
///// Then, it uses the formula $$c(u,v) = \sum_{x\in N(u)} d^{<x}(v)$$
//fn crossing_values(graph: &Graph) -> HashMap<(usize, usize), usize> {
//    
//    let mut d_less_than_x: HashMap<(usize,usize), usize> = HashMap::new();
//
//    // #1 compute all d^<x(u) O(|X||Y|)
//    for u in &graph.bnodes {
//        let mut d: usize = 0;
//        for x in &graph.anodes {
//            d_less_than_x.insert((u.id,x.id),d);
//            if u.neighbors.contains(&x.id) {
//                d += 1;
//            }
//        }
//    }
//
//    // #2 compute all crossing values
//    let mut crossing_dict: HashMap<(usize,usize), usize> = HashMap::new();
//   
//    for u in &graph.bnodes {
//        for v in &graph.bnodes {
//            if u.id==v.id {
//                continue;
//            }
//            for x in &u.neighbors {
//                if let Some(c) = crossing_dict.get(&(u.id,v.id)) {
//                    crossing_dict.insert((u.id,v.id),*c+d_less_than_x.get(&(v.id,*x)).unwrap());
//                }
//                else {
//                    crossing_dict.insert((u.id,v.id),*d_less_than_x.get(&(v.id,*x)).unwrap()); 
//                }
//            }
//        }
//    }
//
//    return crossing_dict;
//}

/// Computes crossing values for orientable pairs only.
/// A pair (u,v) is orientable if neither $r_u\leq l_v$ nor $r_v\leq l_u$.
/// It does so by using the same formula as the function crossing_value,
/// except that $d^{<x}(v)$ is only computed for a $x\in [l_v,r_v]$.
/// Then, c(u,v) is only computed for orientable pairs.
/// If a value of $d^{<x}(v)$ is requested for x smaller than $l_v$ or
/// larger than $r_v$ then it is computed on the flies with its trivial values.
fn orientable_crossing_values(graph: &Graph) -> HashMap<(usize, usize), usize>  {

    let mut d_less_than_x: HashMap<(usize,usize), usize> = HashMap::new();

    // #1 compute d^<x(u) values that matter
    for u in &graph.bnodes {
        let mut d: usize = 0;
        for x in (u.left)..=(u.right) {
            d_less_than_x.insert((u.id,x),d);
            if u.neighbors.contains(&x) {
                d += 1;
            }
        }
    }

    // #1.5 prep work for orientable pairs
    let mut left_end: HashMap<usize,usize> = HashMap::new();
    let mut right_end: HashMap<usize,usize> = HashMap::new();
    let mut neighbors: HashMap<usize, Vec<usize>> = HashMap::new();
    for u in &graph.bnodes {
        left_end.insert(u.id, u.left);
        right_end.insert(u.id, u.right);
        neighbors.insert(u.id, u.neighbors.clone());
    }

    let mut orientable_pairs: HashSet<(usize,usize)> = HashSet::new();
    let mut active_vertices: HashSet<usize> = HashSet::new();
    // #2 computing orientable pairs
    for x in &graph.anodes {
        for u in &x.neighbors {
            if x.id==*left_end.get(u).unwrap() {
                active_vertices.insert(*u);
            }
            if x.id==*right_end.get(u).unwrap() {
                active_vertices.remove(u);
            }
            for w in &active_vertices {
                orientable_pairs.insert((*u, *w));
                orientable_pairs.insert((*w, *u));
            }
        }
    }

    // #3 compute crossing values for orientable pairs
    let mut crossing_dict: HashMap<(usize,usize), usize> = HashMap::new();
   
    for (u,v) in orientable_pairs {
        if u==v {
            continue;
        }
        for x in neighbors.get(&u).unwrap() {

            // reconstructing d<x(v)
            let mut d = 0;
            if x > right_end.get(&v).unwrap() {
                d = neighbors.get(&v).unwrap().len(); 
            }
            else {
                if x > left_end.get(&v).unwrap() {
                    d = *d_less_than_x.get(&(v,*x)).unwrap();
                }
            }

            // adding to crossing value
            if let Some(c) = crossing_dict.get(&(u,v)) {
                crossing_dict.insert((u,v),*c+d);
            }
            else {
                crossing_dict.insert((u,v),d); 
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

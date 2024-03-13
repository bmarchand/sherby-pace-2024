use sherby_pace_2024::*;
use std::path::PathBuf;

#[test]
fn test_twin_pre_processing() {
    let mut path = PathBuf::new();
    path.push("tests/twin.gr");
    let graph: Graph = parse_graph(&path);

    let twin_list = find_twins(&graph);
    println!("{:?}", twin_list);
    assert_eq!(twin_list.len(), 3);

    let graph = process_twins_graph(graph, &twin_list);
    assert_eq!(graph.bnodes.len(), 2);
}

#[test]
fn test_twin_pre_processing2() {
    let mut path = PathBuf::new();
    path.push("tests/twin2.gr");
    let graph: Graph = parse_graph(&path);

    let twin_list = find_twins(&graph);
    println!("{:?}", twin_list);
    assert_eq!(twin_list.len(), 4);

    let graph = process_twins_graph(graph, &twin_list);
    assert_eq!(graph.bnodes.len(), 2);
}

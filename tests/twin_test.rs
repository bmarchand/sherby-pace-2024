use sherby_pace_2024::*;
//use std::path::PathBuf;

#[test]
fn test_twin_pre_processing() {
    let graph: Graph = parse_graph(&"tests/twin.gr".to_string());

    let twin_list = find_twins(&graph);
    println!("{:?}", twin_list);
    assert_eq!(twin_list.len(), 3);

    let graph = process_twins_graph(graph, &twin_list);
    assert_eq!(graph.bnodes.len(), 2);
}

#[test]
fn test_twin_pre_processing2() {
    let graph: Graph = parse_graph(&"tests/twin2.gr".to_string());

    let twin_list = find_twins(&graph);
    println!("{:?}", twin_list);
    assert_eq!(twin_list.len(), 4);

    let graph = process_twins_graph(graph, &twin_list);
    assert_eq!(graph.bnodes.len(), 2);
}

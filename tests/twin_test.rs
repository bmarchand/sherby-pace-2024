//use sherby_pace_2024::*;
//use std::path::PathBuf;
//
//#[test]
//fn test_twin_pre_processing() {
//    let content = std::fs::read_to_string("tests/twin.gr").expect("could not read file");
//    for line in content.lines() {
//        std::io::stdin().read_line(&mut line.to_string()).expect("could not read line");
//    }
//    let graph: Graph = parse_graph();
//
//    let twin_list = find_twins(&graph);
//    println!("{:?}", twin_list);
//    assert_eq!(twin_list.len(), 3);
//
//    let graph = process_twins_graph(graph, &twin_list);
//    assert_eq!(graph.bnodes.len(), 2);
//}
//
//#[test]
//fn test_twin_pre_processing2() {
//    let content = std::fs::read_to_string("tests/twin2.gr").expect("could not read file");
//    for line in content.lines() {
//        std::io::stdin().read_line(&mut line.to_string()).expect("could not read line");
//    }
//    let graph: Graph = parse_graph();
//
//    let twin_list = find_twins(&graph);
//    println!("{:?}", twin_list);
//    assert_eq!(twin_list.len(), 4);
//
//    let graph = process_twins_graph(graph, &twin_list);
//    assert_eq!(graph.bnodes.len(), 2);
//}

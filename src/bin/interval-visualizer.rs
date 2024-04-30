//use clap::Parser;
use sherby_pace_2024::*;
use std::fs::File;
use std::io::Write;

fn main() {
//    let args = Cli::parse();

    let graph_name = std::env::args().nth(1).expect("expecting a path to a graph file.");
    let graph: Graph = parse_graph(&graph_name);
    println!("graph {:?}", graph);

    let filename = graph_name.clone() + ".svg";
    print!("writing to {:?}", filename);
    write_graph_svg(&filename, &graph, true);
}

fn write_graph_svg(filename: &String, graph: &Graph, interval_mode: bool) {
    let unitw: usize = 5;
    let unith: usize = 10;
    let mut maxr: usize = 0;
    let mut ints = Vec::new();

    for node in &graph.bnodes {
        let pl: (usize, usize, usize) = (node.id, node.left as usize, node.right as usize);
        ints.push(pl);
        //print!("{:?}", pl);
    }

    //old syntax I know
    for (_i, _l, r) in &ints {
        if maxr < *r {
            maxr = *r;
        }
    }

    let mut data_file = File::create(filename).expect("Failed to create file");

    write!(data_file, "<!DOCTYPE html>\n").unwrap();
    write!(data_file, "<html>\n").unwrap();
    write!(data_file, "<body>\n").unwrap();

    let mut docheight: usize = 200;

    if interval_mode {
        docheight = ints.len() * unith;
    }

    write!(
        data_file,
        "<svg width='{:?}' height='{:?}' xmlns='http://www.w3.org/2000/svg'>\n",
        maxr * unitw,
        docheight
    )
    .unwrap();

    let mut cpt: usize = 0;

    if interval_mode {
        for (_i, l, r) in &ints {
            write!(data_file,
                "<line x1='{:?}' y1='{:?}' x2='{:?}' y2='{:?}' stroke='black' stroke-width='7' />\n",
                l * unitw, cpt * unith,
                r * unitw + unitw / 2,
                cpt * unith).unwrap();
            cpt += 1;
        }
    } else {
        for node in &graph.bnodes {
            for v in &node.neighbors {
                write!(data_file, "<line x1='{:?}' y1='{:?}' x2='{:?}' y2='{:?}' stroke='black' stroke-width='0.1' />\n", 
    					cpt * unitw, 200, v * unitw, 100).unwrap();
            }
            cpt += 1;
        }
    }

    write!(
        data_file,
        "Sorry, your browser does not support inline SVG.\n"
    )
    .unwrap();
    write!(data_file, "</svg>\n").unwrap();
    write!(data_file, "</body>\n").unwrap();
    write!(data_file, "</html>\n").unwrap();
}

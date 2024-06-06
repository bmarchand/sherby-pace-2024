use assert_cmd::prelude::*;
use std::process::Command;

// all of this file is inspired from https://rust-cli.github.io/book/tutorial/testing.html

#[test]
fn test_complete_4_5() {
    let mut cmd = Command::cargo_bin("sherby-pace-2024").unwrap();
    cmd.args(["< tiny_test_set/complete_4_5.gr","> tiny_test_set/complete_4_5.sol"]);
    cmd.assert().success();

    let sol1 = Command::new("pace2024verifier")
        .args([
            "-c",
            "tiny_test_set/complete_4_5.gr",
            "tiny_test_set/complete_4_5.sol",
        ])
        .output()
        .unwrap();

    let sol2 = Command::new("pace2024verifier")
        .args([
            "-c",
            "tiny_test_set/complete_4_5.gr",
            "tiny_test_set_official_solutions/complete_4_5.sol",
        ])
        .output()
        .unwrap();

    assert_eq!(sol1, sol2);
}

fn test_file(instance: &str, sol1_name: &str, sol2_name: &str) {
    // running our solver. writes solution (order) into sol1.
    let stdin = "< ".to_owned();
    let stdout = "> ".to_owned();
    let mut cmd = Command::cargo_bin("sherby-pace-2024").unwrap();
    cmd.args([stdin+instance,stdout+sol1_name]);
    cmd.assert().success();

    let sol1 = Command::new("pace2024verifier")
        .args(["-c", instance, sol1_name])
        .output()
        .unwrap();

    let sol2 = Command::new("pace2024verifier")
        .args(["-c", instance, sol2_name])
        .output()
        .unwrap();

    assert_eq!(sol1, sol2);
}

fn test_file_no_preprocess(instance: &str, sol1_name: &str, sol2_name: &str) {
    // running our solver. writes solution (order) into sol1.
    let stdin = "< ".to_owned();
    let stdout = "> ".to_owned();
    let mut cmd = Command::cargo_bin("no_preprocess").unwrap();
    cmd.args([stdin+instance,stdout+sol1_name]);
    cmd.assert().success();

    let sol1 = Command::new("pace2024verifier")
        .args(["-c", instance, sol1_name])
        .output()
        .unwrap();

    let sol2 = Command::new("pace2024verifier")
        .args(["-c", instance, sol2_name])
        .output()
        .unwrap();

    assert_eq!(sol1, sol2);
}


#[test]
fn test_cycle_8_shuffled() {
    let mut cmd = Command::cargo_bin("sherby-pace-2024").unwrap();

    cmd.args(["tiny_test_set/cycle_8_shuffled.gr","tiny_test_set/cycle_8_shuffled.sol"]);
    cmd.assert().success();

    let sol1 = Command::new("pace2024verifier")
        .args([
            "-c",
            "tiny_test_set/cycle_8_shuffled.gr",
            "tiny_test_set/cycle_8_shuffled.sol",
        ])
        .output()
        .unwrap();

    let sol2 = Command::new("pace2024verifier")
        .args([
            "-c",
            "tiny_test_set/cycle_8_shuffled.gr",
            "tiny_test_set_official_solutions/cycle_8_shuffled.sol",
        ])
        .output()
        .unwrap();

    assert_eq!(sol1, sol2);
}

#[test]
fn test_low_edges_tiny_test() {
    test_file(
        "tests/motif_test_instances/low_edges_tiny_test.gr",
        "tests/motif_test_instances/low_edges_tiny_test.sol",
        "tests/motif_test_instances/low_edges_tiny_test.true_sol",
        );
}

#[test]
fn test_high_edges_test() {
    test_file(
        "tests/motif_test_instances/high_edges_test.gr",
        "tests/motif_test_instances/high_edges_test.sol",
        "tests/motif_test_instances/high_edges_test.true_sol",
        );
}

#[test]
fn test_high_stars_test() {
    test_file(
        "tests/motif_test_instances/high_stars_test.gr",
        "tests/motif_test_instances/high_stars_test.sol",
        "tests/motif_test_instances/high_stars_test.true_sol",
        );
}

#[test]
fn test_low_edges_test() {
    test_file(
        "tests/motif_test_instances/low_edges_test.gr",
        "tests/motif_test_instances/low_edges_test.sol",
        "tests/motif_test_instances/low_edges_test.true_sol",
        );
}

#[test]
fn test_low_stars_test() {
    test_file(
        "tests/motif_test_instances/low_stars_test.gr",
        "tests/motif_test_instances/low_stars_test.sol",
        "tests/motif_test_instances/low_stars_test.true_sol",
        );
}

#[test]
fn test_low_stars_tiny_test() {
    test_file(
        "tests/motif_test_instances/low_stars_tiny_test.gr",
        "tests/motif_test_instances/low_stars_tiny_test.sol",
        "tests/motif_test_instances/low_stars_tiny_test.true_sol",
        );
}

#[test]
fn test_med_edges_test() {
    test_file(
        "tests/motif_test_instances/med_edges_test.gr",
        "tests/motif_test_instances/med_edges_test.sol",
        "tests/motif_test_instances/med_edges_test.true_sol",
        );
}

#[test]
fn test_med_stars_test() {
    test_file(
        "tests/motif_test_instances/med_stars_test.gr",
        "tests/motif_test_instances/med_stars_test.sol",
        "tests/motif_test_instances/med_stars_test.true_sol",
        );
}

#[test]
fn test_cycle_8_sorted() {
    test_file(
        "tiny_test_set/cycle_8_sorted.gr",
        "tiny_test_set/cycle_8_sorted.sol",
        "tiny_test_set_official_solutions/cycle_8_sorted.sol",
    );
}

#[test]
fn test_grid_9_shuffled() {
    test_file(
        "tiny_test_set/grid_9_shuffled.gr",
        "tiny_test_set/grid_9_shuffled.sol",
        "tiny_test_set_official_solutions/grid_9_shuffled.sol",
    );
}

#[test]
fn test_ladder_4_4_shuffled() {
    test_file(
        "tiny_test_set/ladder_4_4_shuffled.gr",
        "tiny_test_set/ladder_4_4_shuffled.sol",
        "tiny_test_set_official_solutions/ladder_4_4_shuffled.sol",
    );
}

#[test]
fn test_ladder_4_4_sorted() {
    test_file(
        "tiny_test_set/ladder_4_4_sorted.gr",
        "tiny_test_set/ladder_4_4_sorted.sol",
        "tiny_test_set_official_solutions/ladder_4_4_sorted.sol",
    );
}

#[test]
fn test_matching_4_4() {
    test_file(
        "tiny_test_set/matching_4_4.gr",
        "tiny_test_set/matching_4_4.sol",
        "tiny_test_set_official_solutions/matching_4_4.sol",
    );
}

#[test]
fn test_path_9_shuffled() {
    test_file(
        "tiny_test_set/path_9_shuffled.gr",
        "tiny_test_set/path_9_shuffled.sol",
        "tiny_test_set_official_solutions/path_9_shuffled.sol",
    );
}

#[test]
fn test_path_9_sorted() {
    test_file(
        "tiny_test_set/path_9_sorted.gr",
        "tiny_test_set/path_9_sorted.sol",
        "tiny_test_set_official_solutions/path_9_sorted.sol",
    );
}

#[test]
fn test_plane_5_6() {
    test_file(
        "tiny_test_set/plane_5_6.gr",
        "tiny_test_set/plane_5_6.sol",
        "tiny_test_set_official_solutions/plane_5_6.sol",
    );
}

#[test]
fn test_star_6() {
    test_file(
        "tiny_test_set/star_6.gr",
        "tiny_test_set/star_6.sol",
        "tiny_test_set_official_solutions/star_6.sol",
    );
}

#[test]
fn test_tree_6_10() {
    test_file(
        "tiny_test_set/tree_6_10.gr",
        "tiny_test_set/tree_6_10.sol",
        "tiny_test_set_official_solutions/tree_6_10.sol",
    );
}

#[test]
fn test_website_20() {
    test_file(
        "tiny_test_set/website_20.gr",
        "tiny_test_set/website_20.sol",
        "tiny_test_set_official_solutions/website_20.sol",
    );
}

#[test]
fn test_custom_3() {
    test_file(
        "tests/custom3.gr",
        "tests/custom3.sol",
        "tests/custom3_solution.sol",
    );
}

#[test]
fn test_custom_1() {
    test_file(
        "tests/custom1.gr",
        "tests/custom1.sol",
        "tests/custom1_solution.sol",
    );
}

#[test]
fn test_complete_4_5_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/complete_4_5.gr",
        "tiny_test_set/complete_4_5.no_preprocess_sol",
        "tiny_test_set_official_solutions/complete_4_5.sol",
    )
}

#[test]
fn test_cycle_8_shuffled_no_preprocess() {
    let mut cmd = Command::cargo_bin("no_preprocess").unwrap();

    cmd.args(["tiny_test_set/cycle_8_shuffled.gr","tiny_test_set/cycle_8_shuffled.no_preprocess_sol"]);
    cmd.assert().success();

    let sol1 = Command::new("pace2024verifier")
        .args([
            "-c",
            "tiny_test_set/cycle_8_shuffled.gr",
            "tiny_test_set/cycle_8_shuffled.no_preprocess_sol",
        ])
        .output()
        .unwrap();

    let sol2 = Command::new("pace2024verifier")
        .args([
            "-c",
            "tiny_test_set/cycle_8_shuffled.gr",
            "tiny_test_set_official_solutions/cycle_8_shuffled.sol",
        ])
        .output()
        .unwrap();

    assert_eq!(sol1, sol2);
}

#[test]
fn test_cycle_8_sorted_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/cycle_8_sorted.gr",
        "tiny_test_set/cycle_8_sorted.no_preprocess_sol",
        "tiny_test_set_official_solutions/cycle_8_sorted.sol",
    );
}

#[test]
fn test_grid_9_shuffled_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/grid_9_shuffled.gr",
        "tiny_test_set/grid_9_shuffled.no_preprocess_sol",
        "tiny_test_set_official_solutions/grid_9_shuffled.sol",
    );
}

#[test]
fn test_ladder_4_4_shuffled_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/ladder_4_4_shuffled.gr",
        "tiny_test_set/ladder_4_4_shuffled.no_preprocess_sol",
        "tiny_test_set_official_solutions/ladder_4_4_shuffled.sol",
    );
}

#[test]
fn test_ladder_4_4_sorted_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/ladder_4_4_sorted.gr",
        "tiny_test_set/ladder_4_4_sorted.no_preprocess_sol",
        "tiny_test_set_official_solutions/ladder_4_4_sorted.sol",
    );
}

#[test]
fn test_matching_4_4_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/matching_4_4.gr",
        "tiny_test_set/matching_4_4.no_preprocess_sol",
        "tiny_test_set_official_solutions/matching_4_4.sol",
    );
}

#[test]
fn test_path_9_shuffled_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/path_9_shuffled.gr",
        "tiny_test_set/path_9_shuffled.no_preprocess_sol",
        "tiny_test_set_official_solutions/path_9_shuffled.sol",
    );
}

#[test]
fn test_path_9_sorted_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/path_9_sorted.gr",
        "tiny_test_set/path_9_sorted.no_preprocess_sol",
        "tiny_test_set_official_solutions/path_9_sorted.sol",
    );
}

#[test]
fn test_plane_5_6_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/plane_5_6.gr",
        "tiny_test_set/plane_5_6.no_preprocess_sol",
        "tiny_test_set_official_solutions/plane_5_6.sol",
    );
}

#[test]
fn test_star_6_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/star_6.gr",
        "tiny_test_set/star_6.no_preprocess_sol",
        "tiny_test_set_official_solutions/star_6.sol",
    );
}

#[test]
fn test_tree_6_10_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/tree_6_10.gr",
        "tiny_test_set/tree_6_10.no_preprocess_sol",
        "tiny_test_set_official_solutions/tree_6_10.sol",
    );
}

#[test]
fn test_website_20_no_preprocess() {
    test_file_no_preprocess(
        "tiny_test_set/website_20.gr",
        "tiny_test_set/website_20.no_preprocess_sol",
        "tiny_test_set_official_solutions/website_20.sol",
    );
}

#[test]
fn test_custom_3_no_preprocess() {
    test_file_no_preprocess(
        "tests/custom3.gr",
        "tests/custom3.no_preprocess_sol",
        "tests/custom3_solution.sol",
    );
}

#[test]
fn test_custom_1_no_preprocess() {
    test_file_no_preprocess(
        "tests/custom1.gr",
        "tests/custom1.no_preprocess_sol",
        "tests/custom1_solution.sol",
    );
}

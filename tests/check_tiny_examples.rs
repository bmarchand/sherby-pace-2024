use assert_cmd::prelude::*;
use std::process::Command;

// all of this file is inspired from https://rust-cli.github.io/book/tutorial/testing.html

#[test]
fn test_complete_4_5() {
    let mut cmd = Command::cargo_bin("sherby-pace-2024").unwrap();
    cmd.arg("tiny_test_set/complete_4_5.gr");
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
    let mut cmd = Command::cargo_bin("sherby-pace-2024").unwrap();
    cmd.arg(instance);
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

    cmd.arg("tiny_test_set/cycle_8_shuffled.gr");
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
        "tiny_test_set/path_9_shuffled.gr",
        "tiny_test_set/path_9_shuffled.sol",
        "tiny_test_set_official_solutions/path_9_shuffled.sol",
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

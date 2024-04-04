Rust implementation of https://link.springer.com/article/10.1007/s00453-014-9872-x as a starting point for PACE 2024

Resources for learning Rust:  https://www.rust-lang.org/

PACE 2024: https://pacechallenge.org/2024/

Choice of Rust: a short way to describe Rust is that it is similar to C and C++, but with an aggressive, well-designed compiler that forbids you from doing anything dangerous. It is designed to be a replacement of C++. Yet the syntax includes some nice Python-looking stuff. Examples of Rust implementations of a lot of common algorithms: https://github.com/TheAlgorithms/Rust/blob/master/DIRECTORY.md. While it depends on the implementation: Rust is basically as fast as C and C++: https://benchmarksgame-team.pages.debian.net/benchmarksgame/fastest/rust.html

The repository follows the standard Rust structure, as created automatically by Cargo.
The default binary (src/main.rs) solves the problem given a graph. An additional binary (src/bin/interval-visualizer.rs)
yields a visualization of intervals.

Upon cloning the repository, and with Rust installed, running the code on an example should look like:

```
cargo run -- tests/tiny_test_set/cycle_8_shuffled.gr
```

The tests (running on tiny_test_set and cheking that the results are correct) can be run with:
```
cargo test
```

You can use the interval visualizer using:

```
cargo run --bin interval-visualizer tiny_test_set/complete_4_5.gr  
```

You can generate the documentation using 

```
cargo doc --open
```

The documentation is then in target/doc/sherby-pace-2024.

For an actual test of the solvers on the exact data set, you can do:
```
cargo build --release
sudo ./exact_public_exec.sh
```
It compiles with optimization flags, and runs the binary on all instances of the exact
track, with a memory limit of 8GB and a time limit that can be adjusted (in seconds)
in exact_public_exec.sh

For the cutwidth data set, there is a similar script (cutwidth_public_exec.sh).

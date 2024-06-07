Submission for PACE 2024: https://pacechallenge.org/2024/ from the University of Sherbrooke.

To build the solvers:
```
cmake .
make
```

This will trigger the compilation of the SAT-based, C++ part of the solver, based on https://epubs.siam.org/doi/abs/10.1137/1.9781611977561.ch4,
and the rust-part, based on https://link.springer.com/article/10.1007/s00453-014-9872-x.

The repository follows the standard Rust structure, as created automatically by Cargo.
The default binary (src/main.rs) was submitted to the exact track of PACE-2024.
Another binary (src/bin/sherby-cutwidth-hybrid.rs) was submitted to the exact track.
Both first look at the cutwidth of the instance. If it is low enough, the cutwidth-based
solver is called. Otherwise, the SAT solver is called.

After compilation,running the code on an example can be done with:

```
cargo run -- < tests/tiny_test_set/cycle_8_shuffled.gr
```
or 
```
./target/release/sherby-pace-2024 < tests/tiny_test_set/cycle_8_shuffled.gr
```
for the exact track, and
```
./target/release/sherby-cutwidth-hybrid < tests/tiny_test_set/cycle_8_shuffled.gr
```
for the parameterized (cutwidth) track.

The tests (running on tiny_test_set and cheking that the results are correct) can be run with:
```
cargo test
```

Submission for PACE 2024: https://pacechallenge.org/2024/ from the University of Sherbrooke.

**Building the solvers**:
```
cmake .
make
```

This will trigger the compilation of the SAT-based, C++ part of the solver, based on https://epubs.siam.org/doi/abs/10.1137/1.9781611977561.ch4,
and the rust-part, based on https://link.springer.com/article/10.1007/s00453-014-9872-x.

**Dependencies**:
The ``cmake .`` and ``make`` commands above will trigger **cargo**, the standard rust solver. Cargo may require an active internet
connection to fetch automatically some rust dependencies (tempfile, clap, peak_alloc and petgraph) for the project.

**Repository structure and usage**
The repository follows the standard Rust structure, as created automatically by Cargo.
The default binary (src/main.rs) was submitted to the exact track of PACE-2024.
Another binary (src/bin/sherby-cutwidth-hybrid.rs) was submitted to the parameterized track.
Both first look at the cutwidth of the instance. If it is low enough, the cutwidth-based
solver is called. Otherwise, the SAT solver is called.

**Using the solvers:**
```
sherby-exact < exact-public-instances/1.gr
```
for the exact track, and
```
sherby-cutwidth < cutwidth-public/1.gr
```
for the parameterized (cutwidth) track.

The tests (running on tiny_test_set and cheking that the results are correct) can be run with:
```
cargo test
```

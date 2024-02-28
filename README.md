Rust implementation of https://link.springer.com/article/10.1007/s00453-014-9872-x as a starting point for PACE 2024

Resources for learning Rust:  https://www.rust-lang.org/

PACE 2024: https://pacechallenge.org/2024/

Choice of Rust: a short way to describe Rust is that it is similar to C and C++, but with an aggressive, well-designed compiler that forbids you from doing anything dangerous. It is designed to be a replacement of C++. Yet the syntax includes some nice Python-looking stuff. Examples of Rust implementations of a lot of common algorithms: https://github.com/TheAlgorithms/Rust/blob/master/DIRECTORY.md

**Status**
  - implementation of Tamaki-Kobayashi is not yet complete
  - what is an algorithm parameterized by cut-width ? don't know yet.

The repository follows the standard Rust structure, as created automatically by Cargo.

Upon cloning the repository, and with Rust installed, running the code on an example should look like:

```
cargo run -- tiny_test_set/cycle_8_shuffled.gr
```

You can generate the documentation using 

```
cargo doc --open
```

The documentation is then in target/doc/sherby-pace-2024

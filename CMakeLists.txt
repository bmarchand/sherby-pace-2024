cmake_minimum_required(VERSION 3.10)
project(sherby)

add_custom_target(
    sherby-cutwidth ALL
    COMMAND cargo build --release --target=x86_64-unknown-linux-musl
    COMMAND cp target/release/sherby-cutwidth . 
)

add_subdirectory(dfas_v2)

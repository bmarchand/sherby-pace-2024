cmake_minimum_required(VERSION 3.10)
project(sherby)

add_custom_target(
    sherby-cutwidth ALL
    COMMAND cargo build --profile submit --target=x86_64-unknown-linux-musl
    COMMAND cp target/x86_64-unknown-linux-musl/submit/sherby-cutwidth-hybrid sherby-cutwidth 
    COMMAND cp target/x86_64-unknown-linux-musl/submit/sherby-pace-2024 sherby-exact
)

add_subdirectory(dfas_v2)

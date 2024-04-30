cmake_minimum_required(VERSION 3.10)
project(sherby)

add_custom_target(
    sherby-cutwidth ALL
    COMMAND cargo build --release
    COMMAND cp target/release/sherby-cutwidth . 
)

add_subdirectory(dfas_v2)

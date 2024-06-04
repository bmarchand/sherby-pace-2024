cmake_minimum_required(VERSION 3.10)
project(dfas)

set(CMAKE_CXX_STANDARD 20)


#this was a test with clang, but glucose doesn't seem to compile with that
#SET (CMAKE_C_COMPILER:FILEPATH=/usr/bin/clang)
#SET (CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/clang++)



include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if(result)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

set(MYAPP_PGO "")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -Wall -O3 -march=native")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -lz -lm -Wall -O3 -march=native")

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
add_compile_definitions(VERBOSITY=0)
add_subdirectory(EvalMaxSAT/lib/EvalMaxSAT EvalMaxSAT)
#add_subdirectory(EvalMaxSAT/lib/glucose glucose)
#add_subdirectory(EvalMaxSAT/lib/MaLib MaLib)
#add_subdirectory(EvalMaxSAT/lib/cadical cadical)
set(localLibs
        MaLib
        EvalMaxSAT
		glucose
        )
include_directories(EvalMaxSAT/lib/EvalMaxSAT/src)
include_directories(EvalMaxSAT/lib/MaLib/src)
include_directories(EvalMaxSAT/lib/glucose/src)

add_library(dfas_lib SHARED main.cpp)
add_executable(dfas main.cpp)
#target_link_libraries (dfas LINK_PUBLIC glucose MaLib EvalMaxSAT)
target_link_libraries (dfas LINK_PUBLIC glucose MaLib EvalMaxSAT)

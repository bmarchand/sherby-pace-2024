cmake_minimum_required(VERSION 3.16)
project(EvalMaxSAT)

set(CMAKE_CXX_STANDARD 14)

include_directories(lib)
include_directories(lib/EvalMaxSAT/src)
include_directories(lib/glucose/src/glucose/core)
include_directories(lib/glucose/src/glucose/mtl)
include_directories(lib/glucose/src/glucose/simp)
include_directories(lib/glucose/src/glucose/utils)
include_directories(lib/MaLib/src)

add_executable(EvalMaxSAT
        lib/EvalMaxSAT/src/card_oe.cpp
        lib/EvalMaxSAT/src/card_oe.h
        lib/EvalMaxSAT/src/cardincremental.cpp
        lib/EvalMaxSAT/src/cardincremental.h
        lib/EvalMaxSAT/src/EvalMaxSAT.h
        lib/EvalMaxSAT/src/lazyvariable.cpp
        lib/EvalMaxSAT/src/lazyvariable.h
        lib/EvalMaxSAT/src/mcqd.h
        lib/EvalMaxSAT/src/monglucose41.h
        lib/EvalMaxSAT/src/virtualcard.cpp
        lib/EvalMaxSAT/src/virtualcard.h
        lib/EvalMaxSAT/src/virtualmaxsat.h
        lib/EvalMaxSAT/src/virtualsat.cpp
        lib/EvalMaxSAT/src/virtualsat.h
        lib/glucose/src/glucose/core/BoundedQueue.h
        lib/glucose/src/glucose/core/Constants.h
        lib/glucose/src/glucose/core/Dimacs.h
        lib/glucose/src/glucose/core/lcm.cc
        lib/glucose/src/glucose/core/Solver.cc
        lib/glucose/src/glucose/core/Solver.h
        lib/glucose/src/glucose/core/SolverStats.h
        lib/glucose/src/glucose/core/SolverTypes.h
        lib/glucose/src/glucose/mtl/Alg.h
        lib/glucose/src/glucose/mtl/Alloc.h
        lib/glucose/src/glucose/mtl/Clone.h
        lib/glucose/src/glucose/mtl/Heap.h
        lib/glucose/src/glucose/mtl/IntTypes.h
        lib/glucose/src/glucose/mtl/Map.h
        lib/glucose/src/glucose/mtl/Queue.h
        lib/glucose/src/glucose/mtl/Sort.h
        lib/glucose/src/glucose/mtl/Vec.h
        lib/glucose/src/glucose/mtl/VecThreads.h
        lib/glucose/src/glucose/mtl/XAlloc.h
        lib/glucose/src/glucose/simp/Main.cc
        lib/glucose/src/glucose/simp/SimpSolver.cc
        lib/glucose/src/glucose/simp/SimpSolver.h
        lib/glucose/src/glucose/utils/Options.cc
        lib/glucose/src/glucose/utils/Options.h
        lib/glucose/src/glucose/utils/ParseUtils.h
        lib/glucose/src/glucose/utils/System.cc
        lib/glucose/src/glucose/utils/System.h
        lib/MaLib/src/Chrono.h
        lib/MaLib/src/communicationlist.h
        lib/MaLib/src/coutUtil.h
        lib/MaLib/src/main.cpp
        lib/MaLib/src/math.h
        lib/MaLib/src/View.h
        lib/CLI11.hpp
        main.cpp)

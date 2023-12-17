cmake_minimum_required(VERSION 3.27)
project(eikonal_solver)

set(CMAKE_CXX_STANDARD 23)

if(${VERBOSE_SOLVER})
    add_compile_definitions(VERBOSE true)
endif ()
if(${VERBOSE})
    add_compile_definitions(PARSER_VERBOSE true)
    add_compile_definitions(MSHLOADER_VERBOSE true)
endif ()

include_directories("LocalProblem/include")
include_directories(${EIGEN_PATH})
include_directories("include")

file(GLOB SOURCES "src/*.cpp")

set(ENV{OMP_NUM_THREADS} ${OMP_NUM_THREADS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -march=native -funroll-loops -fopenmp")
add_executable(eikonal_solver main.cpp ${SOURCES})

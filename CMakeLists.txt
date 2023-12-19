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

if (NOT DEFINED OMP_NUM_THREADS)
    set(OMP_NUM_THREADS 1)
endif ()

if (NOT DEFINED METHOD)
    set(METHOD FMM)
endif()

if (${METHOD} STREQUAL "FMMP" OR ${METHOD} STREQUAL "fmmp")
    add_compile_definitions(FMMP)
elseif (${METHOD} STREQUAL "FMMP2" OR ${METHOD} STREQUAL "fmmp2")
    add_compile_definitions(FMMP2)
else ()
    add_compile_definitions(FMM)
endif()

set(ENV{OMP_NUM_THREADS} ${OMP_NUM_THREADS})
add_compile_definitions(N_THREADS=${OMP_NUM_THREADS})
add_compile_definitions(${METHOD})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -march=native -funroll-loops -fopenmp")

if (NOT DEFINED ${EXE_NAME})
    set(EXE_NAME "eikonal_solver_${METHOD}")
endif ()

add_executable(${EXE_NAME} main.cpp ${SOURCES})

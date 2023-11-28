cmake_minimum_required(VERSION 3.27)
project(eikonal_test)

set(CMAKE_CXX_STANDARD 23)
include_directories("include")
include_directories("LocalProblem/include")
file(GLOB SOURCES "LocalProblem/*.cpp")
file(GLOB UTILSOURCES "utils/*.cpp")
file(GLOB HEADERS "include/*.h")
file(GLOB HPPEADERS "include/*.hpp")
#also add all cpp file in /
#add all to executable
add_compile_definitions(DIMENSION=3)
#flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -march=native -ffast-math -funroll-loops")
add_executable(eikonal_test main.cpp
        ${SOURCES} ${UTILSOURCES})

target_sources(eikonal_test PRIVATE ${SOURCES} PRIVATE ${UTILSOURCES} PRIVATE ${HEADERS} PRIVATE ${HPPEADERS})

#compile with macro DIMENSION=2
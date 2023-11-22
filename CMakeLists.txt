cmake_minimum_required(VERSION 3.27)
project(eikonel_test)

set(CMAKE_CXX_STANDARD 23)
include_directories("include")
include_directories("LocalProblem/include")
file(GLOB SOURCES "LocalProblem/*.cpp")
#also add all cpp file in /
#add all to executable
add_compile_definitions(DIMENSION=2)
add_executable(eikonel_test main.cpp
        include/PointsEdge.h
        ${SOURCES})

target_sources(eikonel_test PRIVATE ${SOURCES})

#compile with macro DIMENSION=2
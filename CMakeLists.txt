cmake_minimum_required(VERSION 3.27)
project(eikonel_test)

set(CMAKE_CXX_STANDARD 23)
include_directories("include")
add_executable(eikonel_test main.cpp
        include/PointsEdge.h)
file(GLOB SOURCES "LocalProblem/*.cpp")
#also add all cpp file in /
file(GLOB SOURCES ".*cpp")
#add all to executable
target_sources(eikonel_test PRIVATE ${SOURCES})
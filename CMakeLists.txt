cmake_minimum_required(VERSION 3.27)
project(eikonel_test)

set(CMAKE_CXX_STANDARD 23)
include_directories("include")
include_directories("LocalProblem/include")
file(GLOB SOURCES "LocalProblem/*.cpp")
file(GLOB UTILSOURCES "utils/*.cpp")
#also add all cpp file in /
#add all to executable
add_compile_definitions(DIMENSION=3)
add_executable(eikonel_test main.cpp
        include/PointsEdge.h
        ${SOURCES} ${UTILSOURCES})

target_sources(eikonel_test PRIVATE ${SOURCES} PRIVATE ${UTILSOURCES} include/VtkParser.hpp)

#compile with macro DIMENSION=2
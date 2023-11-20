cmake_minimum_required(VERSION 3.27)
project(eikonel_test)

set(CMAKE_CXX_STANDARD 23)

add_executable(eikonel_test main.cpp)
file(GLOB SOURCES "LocalProblem/*.cpp")
#add all to executable
target_sources(eikonel_test PRIVATE ${SOURCES})
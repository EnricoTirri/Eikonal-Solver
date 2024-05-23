cmake_minimum_required(VERSION 3.27)
project(eikonal_solver)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp")

# ------------------- Verbose output Options ------------------ #

if(NOT DEFINED IO_VERBOSE)
    set(IO_VERBOSE "FALSE")
endif ()
if(NOT DEFINED SOLVER_VERBOSE)
    set(SOLVER_VERBOSE "FALSE")
endif ()
if(NOT DEFINED VERBOSE)
    set(VERBOSE "FALSE")
endif ()

if((${IO_VERBOSE} STREQUAL "TRUE") OR (${VERBOSE} STREQUAL "TRUE"))
    add_compile_definitions(PARSER_VERBOSE)
    add_compile_definitions(MSHLOADER_VERBOSE)
endif()

if((${SOLVER_VERBOSE} STREQUAL "TRUE") OR (${VERBOSE} STREQUAL "TRUE"))
    add_compile_definitions(SOLVER_VERBOSE)
endif()

# ------------------------------------------------------------- #

# ----------- Method Selector --------------------------------- #

if (NOT DEFINED METHOD)
    set(METHOD FMM)
endif()

if (${METHOD} STREQUAL "PFMM" OR ${METHOD} STREQUAL "pfmm")
    add_compile_definitions(PFMM) # Enable Patch-Parallelized Fast Marching Method

elseif (${METHOD} STREQUAL "FIMP" OR ${METHOD} STREQUAL "fimp")
    add_compile_definitions(FIMP) # Enable Fast Iterative Method with OpenMP Parallelization

elseif (${METHOD} STREQUAL "FMMP" OR ${METHOD} STREQUAL "fmmp")
    add_compile_definitions(FMMP) # Enable Fast Marching Method with OpenMP Parallelization

elseif (${METHOD} STREQUAL "FIM" OR ${METHOD} STREQUAL "fim")
    add_compile_definitions(FIM) # Enable Fast Iterative Method

else () # DEFAULT
    add_compile_definitions(FMM) # Enable Fast Marching Method

endif()

# ------------------------------------------------------------ #

# ------------- OpenMP Options ------------------------------- #

if (NOT DEFINED OMP_NUM_THREADS)
    set(OMP_NUM_THREADS 1)
endif ()

set(ENV{OMP_NUM_THREADS} ${OMP_NUM_THREADS})
add_compile_definitions(N_THREADS=${OMP_NUM_THREADS})

# ------------------------------------------------------------ #


# --------------- Executable Options ------------------------- #

if (NOT DEFINED ${EXE_NAME})
    set(EXE_NAME "eikonal_solver_${METHOD}")
endif ()

file(GLOB SOURCES "src/*.cpp")
add_executable(${EXE_NAME} main.cpp ${SOURCES})
add_executable("test" test.cpp ${SOURCES})

# ------------------------------------------------------------ #

# ---------- Includes and libraries -------------------------- #

include_directories(include)

## EIGEN
include_directories(${EIGEN_PATH})

## METIS
#/usr/lib/libmetis.so is the path to the metis library need to find a way to make it more general
#set(METIS_LIB "/usr/lib/libmetis.so")
if (NOT DEFINED METIS_LIB)
    set(METIS_LIB "/usr/lib/libmetis.so")
endif ()
target_link_libraries(${EXE_NAME} ${METIS_LIB})

# ----------------------------------------------------------- #



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


# ------------------------------------------------------------ #

# ---------- Includes and libraries -------------------------- #

include_directories(include)

## EIGEN
if (NOT DEFINED EIGEN_PATH)
    set(EIGEN_PATH "/usr/include/eigen3")
endif ()
include_directories(${EIGEN_PATH})

## METIS
#/usr/lib/libmetis.so is the path to the metis library need to find a way to make it more general
#set(METIS_LIB "/usr/lib/libmetis.so")
if (NOT DEFINED METIS_LIB)
    set(METIS_LIB "/usr/lib/libmetis.so")
endif ()

# -------------- Methods definitions -------------------------- #

#first method FMM

add_executable(eikonal_solver_FMM main.cpp ${SOURCES})
target_compile_definitions(eikonal_solver_FMM PRIVATE FMM)
target_link_libraries(eikonal_solver_FMM ${METIS_LIB})


#second method FIMP

add_executable(eikonal_solver_FIMP main.cpp ${SOURCES})
target_compile_definitions(eikonal_solver_FIMP PRIVATE FIMP)
target_link_libraries(eikonal_solver_FIMP ${METIS_LIB})


#third method FMMP

add_executable(eikonal_solver_FMMP main.cpp ${SOURCES})
target_compile_definitions(eikonal_solver_FMMP PRIVATE FMMP)
target_link_libraries(eikonal_solver_FMMP ${METIS_LIB})

#forth method FMMO

add_executable(eikonal_solver_FMMO main.cpp ${SOURCES})
target_compile_definitions(eikonal_solver_FMMO PRIVATE FMMO)
target_link_libraries(eikonal_solver_FMMO ${METIS_LIB})

#fifth method FIM

add_executable(eikonal_solver_FIM main.cpp ${SOURCES})
target_compile_definitions(eikonal_solver_FIM PRIVATE FIM)
target_link_libraries(eikonal_solver_FIM ${METIS_LIB})

#sixth method PFIM

#add_executable(eikonal_solver_PFIM main.cpp ${SOURCES})
#target_compile_definitions(eikonal_solver_PFIM PRIVATE PFIM)
#target_link_libraries(eikonal_solver_PFIM ${METIS_LIB})

#seventh method PFIMC

if(USECUDA)
    enable_language(CUDA)

    include_directories(include-cuda)

    if (NOT DEFINED CUDA_BLOCK_SIZE)
        set(CUDA_BLOCK_SIZE 128)
    endif ()
    add_compile_definitions(BLOCK_SIZE=${CUDA_BLOCK_SIZE})
    #nvcc flags
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_75 -Xcompiler -fopenmp -lcudart  -Xcompiler -O3")
    file(GLOB CUDA_SOURCES "src-cuda/*.c*")
    add_executable(eikonal_solver_PFIMC main.cpp ${SOURCES} ${CUDA_SOURCES})
    target_compile_definitions(eikonal_solver_PFIMC PRIVATE PFIMC)
    target_link_libraries(eikonal_solver_PFIMC ${METIS_LIB})
    set_target_properties(eikonal_solver_PFIMC
            PROPERTIES
            CMAKE_CUDA_STANDARD 23
            CUDA_SEPARABLE_COMPILATION ON
    )
endif ()

# ----------------------------------------------------------- #



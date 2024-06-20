# Eikonal-Negro-Tirri-Santoro
## Overview

This is an implementation of a numerical solver for Eikonal equation.

This is part of a project-work for Advanced Method for Scientific Computing course [@Polimi](https://www.polimi.it/)

### Authors
Project developed by:
- [Negro Giorgio](https://github.com/giorgionegro)
- [Enrico Tirri](https://github.com/EnricoTirri)
- [Santoro Dario](https://github.com/DarioSantoroDS)

### Problem description, Implementation decisions and Result analysis

Can be found into report - `documents/Report.pdf`

### File Structure

* `include` -  headers files where:
    * `Mesh.h` - defines the data structure containing the mesh
    * `VtkParser.hpp` - defines data structure and methods in order to parse a .vtk file
    * `MeshLoader.hpp` - defines methods in order to load a mesh from a file parser
    * `EikonalSolver.hpp` - defines Eikonal global problem solver
    * `EikonalTraits.hpp` - defines some Data Structure for Eikonal problem traits
    * `LocalSolver.hpp` - defines Eikonal local problem solver
    * `OptimizedLocalSolver.hpp` - defines Eikonal Optimized local problem solver
    * `EikonalMath.hpp` - defines common mathematical method implementations

* `include-cuda` - headers files for cuda implemented function wrapping
  * `GlobalSolverKernel.hpp` - defines function wrappers for global solver functions
  * `LocalSolverKernel.hpp` - defines function wrappers for local solver functions

* `src` - source files of headers implementation
* `src-cuda` - source files of cuda function implementation
* `testfiles_torus_tetra` - collection of .vtk file containing tetrahedrical test mesh
* `testfiles_torus_tri` - collection of .vtk file containing triangular test mesh
* `showcase` - collection of reference results images
* `documents` - collection of reference papers and implementation report

### How to build

In order to build the executable, from the root folder run the following commands:

```bash
$ mkdir build
$ cd build
$ cmake .. _FLAGS_
$ make
```
`_FLAGS_` are optional, they can be:
* `-D EIGEN_PATH='/path/to/eigen3'` - specify eigen3 library path
* `-D OMP_NUM_THREADS=_num_threads_` - specify the number of threads used by omp implementation
* `-D VERBOSE=true/false` - specify if you want to enable all verbose options
* `-D SOLVER_VERBOSE=true/false` - specify if show local solver work-messages
* `-D EXE_NAME='_exe_filename_'` - specify the executable filename (default = eikonal_solver_'method_tag')
* ` -D USECUDA=(1 or 0)` - specify if you want to compile the cuda target
* `-D CMAKE_CUDA_COMPILER=` - if you want compile the cuda target it's best to specify the nvcc path, if cuda toolkit is alredy in path it's should not be neccesary
* `-D IO_VERBOSE=true/false` - specify if you want to enable verbose option for vtk file reading and writing
* `-D METIS_LIB='path/to/metis.so'` - specify path to metis shared object
* `-D MAKETEST=(1 or 0)` - specify if compile alternative test targets 

output executables will be:
* `eikonal_solver_FMM` - Fast Marching Method
* `eikonal_solver_FIMP` - Fast Iterative Parallel Method
* `eikonal_solver_FMMO` - Fast Marching Method with Optimized Local solver
* `eikonal_solver_FIM` - Fast Iterative Method
* `eikonal_solver_PFIM` - Patch Fast Iterative Method
* `eikonal_solver_PFIMC` - Patch Fast Iterative Method on GPU

if MAKETEST=1 the output executables will be
* `eikonal_test_FMM` - Fast Marching Method tester
* `eikonal_test_FIMP` - Fast Iterative Parallel Method tester
* `eikonal_test_FMMO` - Fast Marching Method with Optimized Local Solver tester
* `eikonal_test_FIM` - Fast Iterative Method tester
* `eikonal_test_PFIM` - Patch Fast Iterative Method tester
* `eikonal_test_PFIMC` - Patch Fast Iterative Method on GPU tester

### How to run solver

```bash
$ ./eikonal_solver_* input.vtk output.vtk meshdim id1 [id2 ...]
```
where:
* `input.vtk` - is a path/filename to ASCII formatted .vtk input mesh file
* `output.vtk` - is the path/filename where output (mesh + data) will be placed
* `meshdim` - is the mesh size (3=triangular mesh, 4=tetrahedral mesh)
  * **physical dimension is assumed to be always 3,\
    be careful .vtk file to be appropriately formatted**
* `id1` - the .vtk point index of wave starting point [$u(id_1) = 0$]
* `[id2 ...]` - optional list of other starting points

### How to run tester
```bash
$ ./eikonal_test_* input_dir_path output.csv meshdim
```
where:
* `input_dir_path` - is a path to a dir containing **only** ASCII formatted .vtk input mesh files
* `output.csv` - is the path/filename where output (test result) will be placed
* `meshdim` - is the mesh size (3=triangular mesh, 4=tetrahedral mesh)
  * **physical dimension is assumed to be always 3,\
    be careful .vtk file to be appropriately formatted**

### Examples showcase 
##### Triangular Mesh
<img src="showcase/triangles3d.png" width="500"/>

##### Tetrahedral Mesh
<img src="showcase/tetrahedron3d.png" width="500"/>


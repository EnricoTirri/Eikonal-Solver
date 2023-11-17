# Local problem for Eikonal equation #
Based on a modified version of the LineSearch PACS example. It implements a constrained optimization procedure by employing a projected Newton method.

The main driver if the function in `solveEikonalProblem.cpp`. The cpp
marcro *DIMENSION* is used to specify the dimension of the problem (2
or 3). The Makefile will produce the 2D and 3D version of the simple
test contained in main_eikonal.cpp.

A better organization would be to put all instances precompiled in a
library. This is indeed done in the version from which I have
extracted this software, where a python module is created using the
utilities in `pythonBindings/`. The python bindings are here deactivated,
but you can bring them back, that are need pybind11.

Enjoy!

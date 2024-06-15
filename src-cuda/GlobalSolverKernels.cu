#include "GlobalSolverKernels.hpp"
#include "iostream"

namespace Eikonal{

    void allocateAndTransfer(void **dev_ptr, void *host_ptr, unsigned int type_size, unsigned int elem_number){
        cudaMalloc((void **)&(*dev_ptr), type_size * elem_number);
        cudaMemcpy(*dev_ptr, host_ptr, type_size * elem_number, cudaMemcpyHostToDevice);
    }

}

#ifndef EIKONAL_SOLVER_GLOBALSOLVERKERNELS_HPP
#define EIKONAL_SOLVER_GLOBALSOLVERKERNELS_HPP

#include "EikonalTraits.hpp"

namespace Eikonal {

    void allocateAndTransfer(void **dev_ptr, void *host_ptr, unsigned int type_size, unsigned int elem_number);



}

#endif //EIKONAL_SOLVER_GLOBALSOLVERKERNELS_HPP

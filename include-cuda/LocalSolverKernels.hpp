
#ifndef EIKONAL_SOLVER_LOCALSOLVERKERNELS_H
#define EIKONAL_SOLVER_LOCALSOLVERKERNELS_H

#include <array>
#include <EikonalTraits.hpp>

namespace Eikonal {

    template<unsigned int MESH_SIZE>
    using MprimeMatrix = typename TTraits<MESH_SIZE>::MprimeMatrix;

    template<unsigned int MESH_SIZE>
    __device__ double solveLocal(const int &pointref, const MprimeMatrix<MESH_SIZE> &MT, double *valin);
}

#endif //EIKONAL_SOLVER_LOCALSOLVERKERNELS_H

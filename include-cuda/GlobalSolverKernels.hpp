
#ifndef EIKONAL_SOLVER_GLOBALSOLVERKERNELS_HPP
#define EIKONAL_SOLVER_GLOBALSOLVERKERNELS_HPP

#include "EikonalTraits.hpp"

namespace Eikonal {

    template<unsigned int MESH_SIZE>
    using MprimeMatrix = typename TTraits<MESH_SIZE>::MprimeMatrix;

    void allocateAndTransfer(void **dev_ptr, void *host_ptr, unsigned int type_size, unsigned int elem_number);

    template<unsigned int MESH_SIZE>
    void globalSolve(const std::vector<int> &XPatches,
                     const std::vector<int> &patchElementPtr, const std::vector<int> &patchAdjElementIdx,
                     const std::vector<int> &patchNodePtr, const std::vector<int> &patchAdjNodeIdx,
                     const std::vector<int> &elementNodePtr, const std::vector<int> &elementAdjNodeIdx,
                     const std::vector<int> &patchPatchPtr, const std::vector<int> &patchAdjPatchIdx,
                     std::vector<double> &U, const std::vector<MprimeMatrix<MESH_SIZE>> &MprimeList,
                     const std::vector<int> &time_reduction_ptr, int reduction_span,
                     bool *result);
}

#endif //EIKONAL_SOLVER_GLOBALSOLVERKERNELS_HPP

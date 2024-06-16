
#if PFIM

#include "metis.h"
#include "EikonalSolver.hpp"
#include <iostream>

namespace Eikonal {

    template<int MESH_SIZE>
    void EikonalSolver<MESH_SIZE>::print_spec() {
        std::cout << "Patch Fast Iterative Method" << std::endl;
    }

    template<int MESH_SIZE>
    inline void
    partitionMesh(const int &npatch, const Mesh<MESH_SIZE> &mesh, const std::vector<int> &startingPoint,
                  std::vector<int> &adjPatchPtr, std::vector<int> &patchAdjacentPatchList,
                  std::vector<int> &adjPatchElementPointer, std::vector<int> &patchAdjacentElementList,
                  std::vector<int> &startingPatches) {

        // ***************** METIS LIBRARY CALL *********************************************************//
        // Support variables for Metis lib
        auto np = static_cast<idx_t>(npatch);
        auto nn = static_cast<idx_t>(mesh.adjPointPtr.size() - 1);
        auto ne = static_cast<idx_t>(mesh.adjElementPtr.size() - 1);

        std::vector<idx_t> eptr(mesh.adjElementPtr.size());
        std::vector<idx_t> eidx(mesh.elementAdjacentPointList.size());

        std::copy(mesh.adjElementPtr.begin(), mesh.adjElementPtr.end(), eptr.data());
        std::copy(mesh.elementAdjacentPointList.begin(), mesh.elementAdjacentPointList.end(), eidx.data());

        std::vector<idx_t> epart(ne);
        std::vector<idx_t> npart(nn);

        idx_t objval;

        // Call metis lib
        METIS_PartMeshNodal(&ne, &nn, eptr.data(), eidx.data(), nullptr, nullptr, &np, nullptr, nullptr, &objval, epart.data(), npart.data());
        // ********************************************************************************************** //

        // *************** PATCH TO ELEMENT ADJACENT LIST CREATION *********************************//
        // Count elements for patch
        adjPatchElementPointer.resize(npatch + 1);

        for (const idx_t &pnum: epart) {
            adjPatchElementPointer[pnum + 1]++;
        }

        // Create support vector for creation of patch-element adjacent list
        std::vector<idx_t> tempAdjPtr(adjPatchElementPointer.size());
        std::copy(adjPatchElementPointer.begin(), adjPatchElementPointer.end(), tempAdjPtr.data());

        // Create patch-element adjacent pointers list
        for (int i = 0; i < np - 1; ++i) {
            adjPatchElementPointer[i + 1] += adjPatchElementPointer[i];
        }

        // Create patch-element adjacent index list
        patchAdjacentElementList.resize(mesh.adjElementPtr.size() - 1);
        for (int i = 0; i < epart.size(); ++i) {
            patchAdjacentElementList[tempAdjPtr[epart[i]]] = i;
            tempAdjPtr[epart[i]]++;
        }

        // **********************************************************************************//


        // ************** PATCH TO PATCH ADJACENT LIST CREATION ***********************************//

        // vector for cross-checking adjacency
        std::vector<std::vector<bool>> adjacentPatchChecks(npatch);
        for(std::vector<bool> &t : adjacentPatchChecks) {
            t.resize(npatch);
            std::fill(t.begin(), t.end(),false);
        }

        // support vector for count patch-belonging nodes
        std::vector<int> counter(npatch);

        // iterate over all elements
        for (int eleId = 0; eleId < mesh.adjElementPtr.size() - 1; ++eleId) {
            int pStart = mesh.adjElementPtr[eleId];
            int pEnd = mesh.adjElementPtr[eleId + 1];

            std::fill(counter.begin(), counter.end(), 0);

            // count how many points for each patch the element contains
            for (int i = pStart; i < pEnd; ++i)
                counter[npart[mesh.elementAdjacentPointList[i]]]++;

            int refPatch = -1;
            int adjPatch = -1;
            for (int i = 0; i < npatch; ++i) {
                // if element has a face completely in the same patch, then that patch is the "reference" one
                if (counter[i] == (MESH_SIZE - 1)) {
                    refPatch = i;
                }else if(counter[i] > 0){
                    // In case a "reference" patch exists, then only one point will be part of another patch
                    adjPatch = i;
                }
            }

            // sign cross-check if reference exists
            if(refPatch>0) {
                adjacentPatchChecks[refPatch][adjPatch] = true;
                adjacentPatchChecks[adjPatch][refPatch] = true;
            }
        }

        // for each patch count adjacent patches
        adjPatchPtr.resize(npatch+1);
        for(int i=0; i<npatch; ++i){
            int count = 0;
            for(int j=0; j<npatch; ++j)
                if(adjacentPatchChecks[i][j]) count++;
            adjPatchPtr[i+1] = count;
        }

        // create patch to patch adjacent pointer list
        for(int i=0; i<npatch; ++i)
            adjPatchPtr[i+1] += adjPatchPtr[i];

        std::vector<int> tempAdjIdxCount(npatch + 1);
        std::copy(adjPatchPtr.begin(), adjPatchPtr.end(), tempAdjIdxCount.data());

        // Create patch-element adjacent index list
        patchAdjacentPatchList.resize(adjPatchPtr[npatch]);
        for (int i = 0; i < npatch; ++i) {
            for(int j = 0; j < npatch; ++j)
                if(adjacentPatchChecks[i][j]){
                    patchAdjacentPatchList[tempAdjIdxCount[i]] = j;
                    tempAdjIdxCount[i]++;
                }
        }
        // ************************************************************************************** //

        // ************************* Retrieve starting patches ********************************** //

        startingPatches.clear();
        std::vector<bool> tActive(npatch, false);
        for(int i=0; i<startingPoint.size(); ++i){
            int nodeId = startingPoint[i];
            int patchId = npart[nodeId];
            tActive[patchId] = true;
        }

        for(int i=0; i<tActive.size(); ++i)
            if(tActive[i])
                startingPatches.push_back(i);

        // ************************************************************************************** //
    }


    template<int MESH_SIZE>
    inline bool EikonalSolver<MESH_SIZE>::solve(std::vector<double> &U,
                                                const std::vector<int> &X,
                                                const Mesh<MESH_SIZE> &data) {

        std::vector<int> adjPatchPatchPtr;
        std::vector<int> patchAdjPatchIdx;
        std::vector<int> adjPatchElePtr;
        std::vector<int> patchAdjEleIdx;
        std::vector<int> XPatches;

        int nPatches = 5;

        partitionMesh<MESH_SIZE>(nPatches, data, X, adjPatchPatchPtr, patchAdjPatchIdx, adjPatchElePtr, patchAdjEleIdx, XPatches);


        return true;
    }


}


#endif

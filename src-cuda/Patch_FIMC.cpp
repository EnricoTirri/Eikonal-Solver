
#if PFIMC

#include "metis.h"
#include "EikonalSolver.hpp"
#include <iostream>
#include "GlobalSolverKernels.hpp"
#include "OptimizedLocalSolver.hpp"

namespace Eikonal {

    template<int MESH_SIZE>
    void EikonalSolver<MESH_SIZE>::print_spec() {
        std::cout << "Patch Fast Iterative Method - CUDA accelerated" << std::endl;
    }

    template<int MESH_SIZE>
    inline void
    partitionMesh(const int &npatch, const Mesh<MESH_SIZE> &mesh, const std::vector<int> &startingPoint,
                  std::vector<int> &adjPatchPtr, std::vector<int> &patchAdjacentPatchList,
                  std::vector<int> &adjPatchElementPointer, std::vector<int> &patchAdjacentElementList,
                  //std::vector<int> &adjPatchNodePointer, std::vector<int> &patchAdjacentNodeList,
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

//        // *************** PATCH TO NODE ADJACENT LIST CREATION *********************************//
//        // Count elements for patch
//        adjPatchNodePointer.resize(npatch + 1);
//
//        for (const idx_t &pnum: npart) {
//            adjPatchNodePointer[pnum + 1]++;
//        }
//
//        // Create support vector for creation of patch-element adjacent list
//        std::vector<idx_t> tempAdjPtr2(adjPatchNodePointer.size());
//        std::copy(adjPatchNodePointer.begin(), adjPatchNodePointer.end(), tempAdjPtr2.data());
//
//        // Create patch-element adjacent pointers list
//        for (int i = 0; i < np - 1; ++i) {
//            adjPatchNodePointer[i + 1] += adjPatchNodePointer[i];
//        }
//
//        // Create patch-element adjacent index list
//        patchAdjacentNodeList.resize(mesh.adjPointPtr.size() - 1);
//        for (int i = 0; i < npart.size(); ++i) {
//            patchAdjacentNodeList[tempAdjPtr2[npart[i]]] = i;
//            tempAdjPtr2[npart[i]]++;
//        }
//
//        // **********************************************************************************//


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
//        std::vector<int> adjPatchNodePtr;
//        std::vector<int> patchAdjNodeIdx;
        std::vector<int> XPatches;

        int nPatches = 5;

        partitionMesh<MESH_SIZE>(nPatches, data, X, adjPatchPatchPtr, patchAdjPatchIdx, adjPatchElePtr, patchAdjEleIdx,/* adjPatchNodePtr, patchAdjNodeIdx, */XPatches);

        //Initialize host MPrimeMatrices
        std::vector<typename TTraits<MESH_SIZE>::MprimeMatrix> MPrimeMatrixPerPatch(patchAdjEleIdx.size());

        Traits::VelocityM M;
        M << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;

        std::array<Traits::Point, MESH_SIZE> points;
        for(int j=0; j<patchAdjEleIdx.size(); ++j){
            int eleId = patchAdjEleIdx[j];
            int pointRangeStart = data.adjElementPtr[eleId];
            int pointRangeEnd = data.adjElementPtr[eleId + 1];

            for(int k = pointRangeStart; k < pointRangeEnd; ++k){
                points[k-pointRangeStart] = data.points[data.elementAdjacentPointList[pointRangeStart]];
            }

            MPrimeMatrixPerPatch[j] = OptimizedLocalSolver<MESH_SIZE>(0,0,M,points).getMprimeMatrix();
        }


        //Extract max adjElement for node
        int maxRange = 0;
        for(int i=0; i<data.adjPointPtr.size()-1; ++i){
            int range = data.adjPointPtr[i+1] - data.adjPointPtr[i];
            if (range > maxRange) maxRange = range;
        }

        int n_nodes = data.points.size();
        int n_elements = data.adjElementPtr.size() - 1;

        std::vector<int> temp_counters(data.adjPointPtr.size() - 1);

        std::vector<int> time_pointers_ids(maxRange * n_nodes);

        int p=0;
        for(int i=0; i<n_elements; ++i){
            int eleId = patchAdjEleIdx[i];
            int ptRS = data.adjElementPtr[eleId];
            int ptRE = data.adjElementPtr[eleId + 1];

            for(int j=ptRS; j<ptRE; ++j){
                int ptId = data.elementAdjacentPointList[j];
                int count = temp_counters[ptId];
                time_pointers_ids[p] = n_nodes * count + i;
                p++;
                count++;
                temp_counters[ptId] = count;
            }
        }


        //INITIALIZE ALL CUDA DATA STRUCTURES
        double *U_dev;
        int *adjPatchPatchPtr_dev;
        int *patchAdjPatchIdx_dev;
        int *adjPatchElePtr_dev;
        int *patchAdjEleIdx_dev;
        int *XPatches_dev;


        bool globalconverge = true;
        while(globalconverge){




        }




        return true;
    }


}


#endif


#if PFIM

#include "metis.h"
#include "EikonalSolver.hpp"
#include <iostream>
#include "EikonalTraits.hpp"
#include "EikonalHeapComparator.cpp"
// #include "LocalSolver.hpp"
#include "OptimizedLocalSolver.hpp"
#include <chrono>

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
            //   std::vector<int> &adjPatchNodePointer, std::vector<int> &patchAdjacentNodeList,
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
        METIS_PartMeshNodal(&ne, &nn, eptr.data(), eidx.data(), nullptr, nullptr, &np, nullptr, nullptr, &objval,
                            epart.data(), npart.data());

        // ********************************************************************************************** //

        // *************** PATCH TO ELEMENT ADJACENT LIST CREATION *********************************//
        // Count elements for patch
        adjPatchElementPointer.resize(npatch + 1);

        for (const idx_t &pnum: epart) {
            adjPatchElementPointer[pnum + 1]++;
        }

        // Create patch-element adjacent pointers list
        for (int i = 0; i < np; ++i) {
            adjPatchElementPointer[i + 1] += adjPatchElementPointer[i];
        }

        // Create support vector for creation of patch-element adjacent list
        std::vector<idx_t> tempAdjPtr(adjPatchElementPointer.size());
        std::copy(adjPatchElementPointer.begin(), adjPatchElementPointer.end(), tempAdjPtr.data());

        // Create patch-element adjacent index list
        patchAdjacentElementList.resize(mesh.adjElementPtr.size() - 1);
        for (int i = 0; i < epart.size(); ++i) {
            patchAdjacentElementList[tempAdjPtr[epart[i]]] = i;
            tempAdjPtr[epart[i]]++;
        }

        // **********************************************************************************//

        // // *************** PATCH TO NODE ADJACENT LIST CREATION *********************************//
        // // Count elements for patch
        // adjPatchNodePointer.resize(npatch + 1);

        // for (const idx_t &pnum : npart)
        // {
        //     adjPatchNodePointer[pnum + 1]++;
        // }

        // // Create patch-element adjacent pointers list
        // for (int i = 0; i < np; ++i)
        // {
        //     adjPatchNodePointer[i + 1] += adjPatchNodePointer[i];
        // }

        // // Create support vector for creation of patch-element adjacent list
        // std::vector<idx_t> tempAdjPtr2(adjPatchNodePointer.size());
        // std::copy(adjPatchNodePointer.begin(), adjPatchNodePointer.end(), tempAdjPtr2.data());

        // // Create patch-element adjacent index list
        // patchAdjacentNodeList.resize(mesh.adjPointPtr.size() - 1);
        // for (int i = 0; i < npart.size(); ++i)
        // {
        //     patchAdjacentNodeList[tempAdjPtr2[npart[i]]] = i;
        //     tempAdjPtr2[npart[i]]++;
        // }

        // **********************************************************************************//

        // ************** PATCH TO PATCH ADJACENT LIST CREATION ***********************************//

        // vector for cross-checking adjacency
        std::vector<std::vector<bool>> adjacentPatchChecks(npatch);
        for (std::vector<bool> &t: adjacentPatchChecks) {
            t.resize(npatch);
            std::fill(t.begin(), t.end(), false);
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
                } else if (counter[i] > 0) {
                    // In case a "reference" patch exists, then only one point will be part of another patch
                    adjPatch = i;
                }
            }

            // sign cross-check if reference exists
            if (refPatch > 0) {
                adjacentPatchChecks[refPatch][adjPatch] = true;
                adjacentPatchChecks[adjPatch][refPatch] = true;
            }
        }

        // for each patch count adjacent patches
        adjPatchPtr.resize(npatch + 1);
        for (int i = 0; i < npatch; ++i) {
            int count = 0;
            for (int j = 0; j < npatch; ++j)
                if (adjacentPatchChecks[i][j])
                    count++;
            adjPatchPtr[i + 1] = count;
        }

        // create patch to patch adjacent pointer list
        for (int i = 0; i < npatch; ++i)
            adjPatchPtr[i + 1] += adjPatchPtr[i];

        std::vector<int> tempAdjIdxCount(npatch + 1);
        std::copy(adjPatchPtr.begin(), adjPatchPtr.end(), tempAdjIdxCount.data());

        // Create patch-element adjacent index list
        patchAdjacentPatchList.resize(adjPatchPtr[npatch]);
        for (int i = 0; i < npatch; ++i) {
            for (int j = 0; j < npatch; ++j)
                if (adjacentPatchChecks[i][j]) {
                    patchAdjacentPatchList[tempAdjIdxCount[i]] = j;
                    tempAdjIdxCount[i]++;
                }
        }
        // ************************************************************************************** //

        // ************************* Retrieve starting patches ********************************** //

        startingPatches.clear();
        std::vector<bool> tActive(npatch, false);
        for (int i = 0; i < startingPoint.size(); ++i) {
            int nodeId = startingPoint[i];
            int patchId = npart[nodeId];
            tActive[patchId] = true;
        }

        for (int i = 0; i < tActive.size(); ++i)
            if (tActive[i])
                startingPatches.push_back(i);

        // ************************************************************************************** //
    }

    template<int MESH_SIZE>
    bool EikonalSolver<MESH_SIZE>::solve(std::vector<double> &U,
                                         const std::vector<int> &X,
                                         const Mesh<MESH_SIZE> &data) {

        auto start = std::chrono::high_resolution_clock::now();

        // prepare partitioning mesh
        std::vector<int> adjPatchPatchPtr;
        std::vector<int> patchAdjPatchIdx;
        std::vector<int> adjPatchElePtr;
        std::vector<int> patchAdjEleIdx;
        std::vector<int> adjPatchNodePtr;
        std::vector<int> patchAdjNodeIdx;
        std::vector<int> XPatches;

        int patchSize = 512;
        int n_elements = data.adjElementPtr.size() - 1;

        int nPatches = (n_elements + patchSize) / patchSize;

        nPatches += (nPatches == 1) ? 1 : 0;

        // partition mesh
        partitionMesh<MESH_SIZE>(nPatches, data, X, adjPatchPatchPtr, patchAdjPatchIdx, adjPatchElePtr,
                                 patchAdjEleIdx, /*adjPatchNodePtr, patchAdjNodeIdx,*/ XPatches);

        // prepare patch convergence vector
        std::vector<bool> convergence(nPatches, false);

        // prepare solution vector
        U.resize(data.points.size());
        std::fill(U.begin(), U.end(), MAXF);
        for (const int &pntID: X) {
            U[pntID] = 0.0;
        }

        // precompute all Optimized Local Solver
        std::vector<OptimizedLocalSolver<MESH_SIZE>> solvers;
        Traits::VelocityM M;
        M << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;
        std::array<Traits::Point, MESH_SIZE> points;
        for (const int &i: patchAdjEleIdx) {
            int ptrStart = data.adjElementPtr[i];
            int ptrEnd = data.adjElementPtr[i + 1];
            for (int j = 0; j < MESH_SIZE; ++j) {
                points[j] = data.points[data.elementAdjacentPointList[j + ptrStart]];
            }
            solvers.push_back(OptimizedLocalSolver<MESH_SIZE>{1, 10e-6, M, points});
        }

        // init active list
        std::vector<bool> activePatchesList(nPatches, false);
        for (const int &i: XPatches) {
            activePatchesList[i] = true;
        }
        bool isActivePatchListEmpty = false;

        auto end = std::chrono::high_resolution_clock::now();
        this->prepare = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        start = std::chrono::high_resolution_clock::now();

        // number of iterations for each convergence cycle
        const int iterations = 2;

#pragma omp parallel
        {
            while (!isActivePatchListEmpty) {

                // init convergence vector
#pragma omp single
                {
                    convergence.assign(convergence.size(), false);
                }
#pragma omp barrier

                // Update all active patches
#pragma omp for
                for (int activePatch = 0; activePatch < nPatches; activePatch++) {

                    if (!activePatchesList[activePatch])
                        continue;

                    for (int it = 0; it < iterations; it++) {

                        // iterate on all elements of the patch
                        for (int activeTriangleIdx = 0; activeTriangleIdx < adjPatchElePtr[activePatch + 1] -
                                                                            adjPatchElePtr[activePatch]; activeTriangleIdx++){

                            int eleId = patchAdjEleIdx[activeTriangleIdx + adjPatchElePtr[activePatch]];

                            int pointRangeStart = data.adjElementPtr[eleId];
                            int pointRangeEnd = data.adjElementPtr[eleId + 1];

                            // retrieve local solver data
                            OptimizedLocalSolver<MESH_SIZE> solver = solvers[eleId];
                            std::array<double, MESH_SIZE> values;
                            for (std::size_t j = 0; j < MESH_SIZE; j++) {
                                int tpId = data.elementAdjacentPointList[pointRangeStart + j];
                                values[j] = U[tpId];
                            }

                            // iterate on all points of elements
                            for (int pi = pointRangeStart; pi < pointRangeEnd; ++pi) {
                                int pointID = data.elementAdjacentPointList[pi];

                                int k = pi - pointRangeStart;

                                double sol = solver(k, values);

                                if (sol >= U[pointID])
                                    continue;

                                values[k] = sol;

#pragma omp atomic write
                                U[pointID] = sol;

                                // if a points has converged then all patch needs to be consideres as converged
                                convergence[activePatch] = 1;
                            }
                        }
                    }
                }

                // Activate neighbour patches of converged patches
#pragma omp for
                for (int activePatch = 0; activePatch < nPatches; activePatch++) {
                    if (!activePatchesList[activePatch])
                        continue;
                    if (convergence[activePatch]) {
                        for (int adjPatchdix = 0; adjPatchdix < adjPatchPatchPtr[activePatch + 1] -
                                                                adjPatchPatchPtr[activePatch]; adjPatchdix++) {
                            activePatchesList[patchAdjPatchIdx[adjPatchPatchPtr[activePatch] + adjPatchdix]] = true;
                        }
                    }
                }

                // Do a single iteration to check if activated patches are going to really converge
#pragma omp for
                for (int activePatch = 0; activePatch < nPatches; activePatch++) {
                    if (!activePatchesList[activePatch])
                        continue;

                    for (int iterations = 0; iterations < 1; iterations++) {

                        for (int activeTriangleIdx = 0; activeTriangleIdx < adjPatchElePtr[activePatch + 1] -
                                                                            adjPatchElePtr[activePatch]; activeTriangleIdx++) {
                            int eleId = patchAdjEleIdx[activeTriangleIdx + adjPatchElePtr[activePatch]];
                            int pointRangeStart = data.adjElementPtr[eleId];
                            int pointRangeEnd = data.adjElementPtr[eleId + 1];

                            OptimizedLocalSolver<MESH_SIZE> solver = solvers[eleId];

                            std::array<double, MESH_SIZE> values;
                            for (std::size_t j = 0; j < MESH_SIZE; j++) {
                                int tpId = data.elementAdjacentPointList[pointRangeStart + j];
                                values[j] = U[tpId];
                            }

                            for (int pi = pointRangeStart; pi < pointRangeEnd; ++pi) {
                                int pointID = data.elementAdjacentPointList[pi];

                                int k = pi - pointRangeStart;

                                double sol = solver(k, values);

                                if (sol >= U[pointID])
                                    continue;
                                values[k] = sol;

#pragma omp atomic write
                                U[pointID] = sol;
                                convergence[activePatch] = 1;
                            }
                        }
                    }
                }

                // re-init activepatchlist
#pragma omp single
                {
                    activePatchesList.assign(activePatchesList.size(), false);
                }

                // re-init isActivePatchListEmpty
#pragma omp single
                {
                    isActivePatchListEmpty = true;
                }

#pragma omp barrier

                // activate patches if converged
#pragma omp for
                for (int i = 0; i < nPatches; i++) {
                    if (convergence[i]) {
                        activePatchesList[i] = true;
                        isActivePatchListEmpty = false;
                    }
                }
            }
        }

        end = std::chrono::high_resolution_clock::now();
        this->compute = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        return true;
    }
}

#endif

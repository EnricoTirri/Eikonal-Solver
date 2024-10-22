#include "GlobalSolverKernels.hpp"
#include "LocalSolverKernels.hpp"
#include "iostream"
#include <thrust/device_ptr.h>
#include <thrust/scan.h>

namespace Eikonal {

#define MAXF 900000.0;

    void allocateAndTransfer(void **dev_ptr, void *host_ptr, unsigned int type_size, unsigned int elem_number) {
        cudaMalloc((void **) &(*dev_ptr), type_size * elem_number);
        cudaMemcpy(*dev_ptr, host_ptr, type_size * elem_number, cudaMemcpyHostToDevice);
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) throw;
    }

    __global__ void
    initTimeReductionList(double *time_reduction_list, int list_size) {
        int blockSize = blockDim.x * blockDim.y;
        int blockId = blockIdx.y * gridDim.x + blockIdx.x;
        int threadId = threadIdx.y * blockDim.x + threadIdx.x;

        int id = blockId * blockSize + threadId;

        if (id < list_size) {
            time_reduction_list[id] = MAXF;
        }
    }

    __global__ void
    activateNeighbours(int *patchPatchPtr, int *patchAdjPatchIdx, int *converged_patch_list,
                       int *converged_patch_list_new, int n_patches) {

        int blockSize = blockDim.x * blockDim.y;
        int blockId = blockIdx.y * gridDim.x + blockIdx.x;
        int threadId = threadIdx.y * blockDim.x + threadIdx.x;
        int id = blockId * blockSize + threadId;

        if (id < n_patches) {
            if (converged_patch_list[id]) {
                int pRS = patchPatchPtr[id];
                int pRE = patchPatchPtr[id + 1];

                for (int i = pRS; i < pRE; ++i) {
                    converged_patch_list_new[patchAdjPatchIdx[i]] = 1;
                }
            }
        }
    }

    __global__ void
    updatePatchValues3(double *U, MprimeMatrix<3> *MprimeList, const int iterations,
                       const int *activePatchList, const int *activePatchList_size,
                       const int *patchElementPtr, const int *patchAdjElementIdx,
                       const int *patchNodePtr, const int *patchAdjNodeIdx,
                       const int *elementNodePtr, const int *elementAdjNodeIdx,
                       double *time_reduction_list, const int *time_reduction_ptr,
                       const int reductionSpan, int *converged_patch_list) {

        // GET BLOCK ID == ACTIVE PATCH POS
        int blockId = blockIdx.y * gridDim.x + blockIdx.x;

        // CHECK IF BLOCK ID IS IN RANGE OF ACTIVE PATCHES
        if (blockId < *activePatchList_size) {
            // GET ACTIVE PATCH ID
            int patchId = activePatchList[blockId];

            // GET THREAD ID == ELEMENT POS IN PATCH
            int threadId = threadIdx.y * blockDim.x + threadIdx.x;

            // GET NUMBER OF NODE AND ELEMENTS OF THE PATCH
            int patchEleStart = patchElementPtr[patchId];
            int patchEleEnd = patchElementPtr[patchId + 1];
            int patchEleSize = patchEleEnd - patchEleStart;

            int patchNodeStart = patchNodePtr[patchId];
            int patchNodeEnd = patchNodePtr[patchId + 1];
            int patchNodeSize = patchNodeEnd - patchNodeStart;

            // RETRIEVE MAX BETWEEN N_ELEMENT AND N_NODES (needed to dispatch threads among elements and nodes)
            int maxn = patchEleSize > patchNodeSize ? patchEleSize : patchNodeSize;

            // IF THREAD IS IN RANGE OF MAX
            if (threadId < maxn) {

                int elePos;
                MprimeMatrix<3> MT;
                int eleId;
                int eleStart;
                int eleEnd;
                double times[3];

                // IF THREAD IS IN RANGE OF ELEMENTS
                if (threadId < patchEleSize) {
                    // get element position in per patch element sorted list
                    elePos = patchEleStart + threadId;

                    // get MprimeMatrix associated to element
                    MT = MprimeList[elePos];

                    // get element Id and range of points
                    eleId = patchAdjElementIdx[elePos];
                    eleStart = elementNodePtr[eleId];
                    eleEnd = elementNodePtr[eleId + 1];
                }

                int nodePos;
                int nodeId;

                // IF THREAD IS IN RANGE OF NODES
                if (threadId < patchNodeSize) {
                    // get node position in per path node sorted list
                    nodePos = patchNodeStart + threadId;
                    // get node id
                    nodeId = patchAdjNodeIdx[nodePos];
                    // init node as not converged
                    //converged_reduction_list[nodePos] = 0;
                }

                if (threadId == 0)
                    converged_patch_list[patchId] = 0;

                // MAIN ITERATION
                for (int iteration = 0; iteration < iterations; ++iteration) {

                    __syncthreads();

                    // IF THREAD IN RANGE OF ELEMENTS
                    if (threadId < patchEleSize) {

                        // RETRIEVE CURRENT TIMES OF VERTICES
                        for (int i = eleStart; i < eleEnd; ++i) {
                            // retrieve current time associated to each point
                            int pointId = elementAdjNodeIdx[i];
                            times[i - eleStart] = U[pointId];
                        }

                        // UPDATE ALL TIMES OF VERTICES
                        for (int i = 0; i < 3; ++i) {

                            // solve for each vertex of element
                            double sol = solveLocal<3>(i, MT, times);

                            // if locally converged
                            if (times[i] > sol) {
                                // locally assign solution
                                times[i] = sol;
                                // place result in the reduction list for later reconciliation of values
                                time_reduction_list[time_reduction_ptr[elePos * 3 + i]] = sol;
                            }
                        }
                    }


                    // WAIT FOR ALL THREADS TO HAVE CALCULATED TIMES OF ELEMENTS
                    __syncthreads();


                    // IF THREAD IN RANGE OF NODES => REDUCTION OF TIME ON NODES OF PATCH
                    if (threadId < patchNodeSize) {

                        // init reduction variable on first of reduction node values
                        int listPos = nodeId * reductionSpan;

                        double reducedValue = time_reduction_list[listPos];
                        // over other reduction node values
                        for (int i = 1; i < reductionSpan; ++i)
                            // apply reduction if i-th value is smaller
                            if (reducedValue > time_reduction_list[listPos + i])
                                reducedValue = time_reduction_list[listPos + i];


                        // check if node has converged
                        if (U[nodeId] > reducedValue) {
                            // set as converged and update node time value
                            converged_patch_list[patchId] = 1;
                            U[nodeId] = reducedValue;
                        }
                    }
                }
            }
        }
    }

    __global__ void
    updatePatchValues4(double *U, MprimeMatrix<4> *MprimeList, const int iterations,
                       const int *activePatchList, const int *activePatchList_size,
                       const int *patchElementPtr, const int *patchAdjElementIdx,
                       const int *patchNodePtr, const int *patchAdjNodeIdx,
                       const int *elementNodePtr, const int *elementAdjNodeIdx,
                       double *time_reduction_list, const int *time_reduction_ptr,
                       const int reductionSpan, int *converged_patch_list) {

        // GET BLOCK ID == ACTIVE PATCH POS
        int blockId = blockIdx.y * gridDim.x + blockIdx.x;

        // CHECK IF BLOCK ID IS IN RANGE OF ACTIVE PATCHES
        if (blockId < *activePatchList_size) {
            // GET ACTIVE PATCH ID
            int patchId = activePatchList[blockId];

            // GET THREAD ID == ELEMENT POS IN PATCH
            int threadId = threadIdx.y * blockDim.x + threadIdx.x;

            // GET NUMBER OF NODE AND ELEMENTS OF THE PATCH
            int patchEleStart = patchElementPtr[patchId];
            int patchEleEnd = patchElementPtr[patchId + 1];
            int patchEleSize = patchEleEnd - patchEleStart;

            int patchNodeStart = patchNodePtr[patchId];
            int patchNodeEnd = patchNodePtr[patchId + 1];
            int patchNodeSize = patchNodeEnd - patchNodeStart;

            // RETRIEVE MAX BETWEEN N_ELEMENT AND N_NODES (needed to dispatch threads among elements and nodes)
            int maxn = patchEleSize > patchNodeSize ? patchEleSize : patchNodeSize;

            // IF THREAD IS IN RANGE OF MAX
            if (threadId < maxn) {

                int elePos;
                MprimeMatrix<4> MT;
                int eleId;
                int eleStart;
                int eleEnd;
                double times[4];

                // IF THREAD IS IN RANGE OF ELEMENTS
                if (threadId < patchEleSize) {
                    // get element position in per patch element sorted list
                    elePos = patchEleStart + threadId;

                    // get MprimeMatrix associated to element
                    MT = MprimeList[elePos];

                    // get element Id and range of points
                    eleId = patchAdjElementIdx[elePos];
                    eleStart = elementNodePtr[eleId];
                    eleEnd = elementNodePtr[eleId + 1];
                }

                int nodePos;
                int nodeId;

                // IF THREAD IS IN RANGE OF NODES
                if (threadId < patchNodeSize) {
                    // get node position in per path node sorted list
                    nodePos = patchNodeStart + threadId;
                    // get node id
                    nodeId = patchAdjNodeIdx[nodePos];
                    // init node as not converged
                    //converged_reduction_list[nodePos] = 0;
                }

                if (threadId == 0)
                    converged_patch_list[patchId] = 0;

                // MAIN ITERATION
                for (int iteration = 0; iteration < iterations; ++iteration) {

                    __syncthreads();

                    // IF THREAD IN RANGE OF ELEMENTS
                    if (threadId < patchEleSize) {

                        // RETRIEVE CURRENT TIMES OF VERTICES
                        for (int i = eleStart; i < eleEnd; ++i) {
                            // retrieve current time associated to each point
                            int pointId = elementAdjNodeIdx[i];
                            times[i - eleStart] = U[pointId];
                        }

                        // UPDATE ALL TIMES OF VERTICES
                        for (int i = 0; i < 4; ++i) {

                            // solve for each vertex of element
                            double sol = solveLocal<4>(i, MT, times);

                            // if locally converged
                            if (times[i] > sol) {
                                // locally assign solution
                                times[i] = sol;
                                // place result in the reduction list for later reconciliation of values
                                time_reduction_list[time_reduction_ptr[elePos * 4 + i]] = sol;
                            }
                        }
                    }


                    // WAIT FOR ALL THREADS TO HAVE CALCULATED TIMES OF ELEMENTS
                    __syncthreads();


                    // IF THREAD IN RANGE OF NODES => REDUCTION OF TIME ON NODES OF PATCH
                    if (threadId < patchNodeSize) {

                        // init reduction variable on first of reduction node values
                        int listPos = nodeId * reductionSpan;

                        double reducedValue = time_reduction_list[listPos];
                        // over other reduction node values
                        for (int i = 1; i < reductionSpan; ++i)
                            // apply reduction if i-th value is smaller
                            if (reducedValue > time_reduction_list[listPos + i])
                                reducedValue = time_reduction_list[listPos + i];


                        // check if node has converged
                        if (U[nodeId] > reducedValue) {
                            // set as converged and update node time value
                            converged_patch_list[patchId] = 1;
                            U[nodeId] = reducedValue;
                        }
                    }
                }
            }
        }
    }

    void checkError(const std::string &message = "") {
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cout << cudaGetErrorName(err) << " " << message << std::endl;
            throw;
        }
    }


    __global__ void
    writeActiveIndices(const int *convergedPatchListNew, const int *scannedList, int *activePatchList, int patches) {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < patches && convergedPatchListNew[idx] > 0) {
            activePatchList[scannedList[idx]] = idx;
        }
    }

#define checkCudaError(err) {\
    if (err != cudaSuccess) {\
        std::cerr << "CUDA Error: " << cudaGetErrorString(err) << std::endl;\
        exit(EXIT_FAILURE);\
    }\
}

    int reduceActiveList(const int *d_convergedPatchListNew, int *d_activePatchList, int *d_activePatchListSize,
                         int patches) {
        // Allocate device memory for the scanned list
        int *d_scannedList;
        checkCudaError(cudaMalloc(&d_scannedList, patches * sizeof(int)));

        // Use Thrust device pointers for Thrust operations
        thrust::device_ptr<const int> d_convergedPatchListNew_ptr(d_convergedPatchListNew);
        thrust::device_ptr<int> d_scannedList_ptr(d_scannedList);
        // Perform an exclusive scan on the convergedPatchListNew array directly
        thrust::exclusive_scan(d_convergedPatchListNew_ptr, d_convergedPatchListNew_ptr + patches, d_scannedList_ptr);


        int *scannedList = (int *) malloc(1 * sizeof(int));
        thrust::copy(d_scannedList_ptr + patches - 1, d_scannedList_ptr + patches, scannedList);

        int activeListSize = scannedList[0];
        //print scanned list for debugging
        int *d_convergedPatchListNew_last = (int *) malloc(1 * sizeof(int));
        thrust::copy(d_convergedPatchListNew_ptr + patches - 1, d_convergedPatchListNew_ptr + patches,
                     d_convergedPatchListNew_last);

        activeListSize += (d_convergedPatchListNew_last[0] > 0 ? 1 : 0);


        checkCudaError(cudaMemcpy(d_activePatchListSize, &activeListSize, sizeof(int), cudaMemcpyHostToDevice));

        // Launch kernel to write active indices
        int threadsPerBlock = 256;
        int blocksPerGrid = (patches + threadsPerBlock - 1) / threadsPerBlock;


        writeActiveIndices<<<blocksPerGrid, threadsPerBlock>>>(d_convergedPatchListNew, d_scannedList,
                                                               d_activePatchList, patches);

        checkCudaError(cudaFree(d_scannedList));
        return activeListSize;
    }

#undef checkCudaError

    template<>
    void globalSolve<3>(const std::vector<int> &XPatches,
                        const std::vector<int> &patchElementPtr, const std::vector<int> &patchAdjElementIdx,
                        const std::vector<int> &patchNodePtr, const std::vector<int> &patchAdjNodeIdx,
                        const std::vector<int> &elementNodePtr, const std::vector<int> &elementAdjNodeIdx,
                        const std::vector<int> &patchPatchPtr, const std::vector<int> &patchAdjPatchIdx,
                        std::vector<double> &U, const std::vector<MprimeMatrix<3>> &MprimeList,
                        const std::vector<int> &time_reduction_ptr, const int reduction_span,
                        bool *result) {

        int block_size = BLOCK_SIZE;

        int n_patches = patchElementPtr.size() - 1;
        int n_elements = patchAdjElementIdx.size();
        int n_nodes = patchAdjNodeIdx.size();

        // INIT CUDA PATCH TO ELEMENT ADJACENT LIST
        int *patchElementPtr_dev, *patchAdjElementIdx_dev;
        {
            allocateAndTransfer((void **) &patchElementPtr_dev, (void *) patchElementPtr.data(),
                                sizeof(int), patchElementPtr.size());
            allocateAndTransfer((void **) &patchAdjElementIdx_dev, (void *) patchAdjElementIdx.data(),
                                sizeof(int), patchAdjElementIdx.size());
        }

        // INIT CUDA PATCH TO NODE ADJACENT LIST
        int *patchNodePtr_dev, *patchAdjNodeIdx_dev;
        {
            allocateAndTransfer((void **) &patchNodePtr_dev, (void *) patchNodePtr.data(),
                                sizeof(int), patchNodePtr.size());
            cudaMemcpy((void *) patchNodePtr.data(), patchNodePtr_dev, sizeof(int) * patchNodePtr.size(),
                       cudaMemcpyDeviceToHost);

            allocateAndTransfer((void **) &patchAdjNodeIdx_dev, (void *) patchAdjNodeIdx.data(),
                                sizeof(int), patchAdjNodeIdx.size());
        }

        // INIT CUDA ELEMENT TO NODE ADJACENT LIST
        int *elementNodePtr_dev, *elementAdjNodeIdx_dev;
        {
            allocateAndTransfer((void **) &elementNodePtr_dev, (void *) elementNodePtr.data(),
                                sizeof(int), elementNodePtr.size());
            allocateAndTransfer((void **) &elementAdjNodeIdx_dev, (void *) elementAdjNodeIdx.data(),
                                sizeof(int), elementAdjNodeIdx.size());
        }

        // INIT CUDA PATCH TO PATCH ADJACENT LIST
        int *patchPatchPtr_dev, *patchAdjPatchIdx_dev;
        {
            allocateAndTransfer((void **) &patchPatchPtr_dev, (void *) patchPatchPtr.data(),
                                sizeof(int), patchPatchPtr.size());
            allocateAndTransfer((void **) &patchAdjPatchIdx_dev, (void *) patchAdjPatchIdx.data(),
                                sizeof(int), patchAdjPatchIdx.size());
        }

        // INIT CUDA MPRIMEMATRIX PER ELEMENT LIST
        MprimeMatrix<3> *MprimeList_dev;
        {
            allocateAndTransfer((void **) &MprimeList_dev, (void *) MprimeList.data(),
                                sizeof(MprimeMatrix<3>), MprimeList.size());
        }

        // INIT CUDA TIME REDUCTION LIST AND POINTER LIST
        int *time_reduction_ptr_dev;
        double *time_reduction_list_dev;
        {
            allocateAndTransfer((void **) &time_reduction_ptr_dev, (void *) time_reduction_ptr.data(),
                                sizeof(int), time_reduction_ptr.size());

            // init all time reduction list at MAXF
            int list_size = reduction_span * n_nodes;
            cudaMalloc((void **) &time_reduction_list_dev, sizeof(double) * list_size);
            int n_blocks = (list_size + block_size) / block_size;
            initTimeReductionList<<<n_blocks, block_size>>>(time_reduction_list_dev, list_size);
        }

        // INIT CUDA SUPPORT LISTS FOR CONVERGED NODES AND PATCHES
        int *converged_patch_list_dev, *converged_patch_list_new_dev;
        {
            cudaMalloc((void **) &converged_patch_list_dev, sizeof(int) * n_patches);
            cudaMalloc((void **) &converged_patch_list_new_dev, sizeof(int) * n_patches);
        }

        // INIT CUDA TIME
        double *U_dev;
        {
            allocateAndTransfer((void **) &U_dev, (void *) U.data(),
                                sizeof(double), U.size());
        }

        // INIT CUDA ACTIVE LIST AND SIZE
        std::vector<int> activePatchList(n_patches);
        std::copy(XPatches.begin(), XPatches.end(), activePatchList.data());

        int activeListSize = XPatches.size();
        int *activePatchList_dev, *activePatchListSize_dev;
        {
            allocateAndTransfer((void **) &activePatchList_dev, (void *) activePatchList.data(),
                                sizeof(int), n_patches);
            allocateAndTransfer((void **) &activePatchListSize_dev, (void *) &activeListSize,
                                sizeof(int), 1);
        }


        std::vector<int> convergedPatchList(n_patches);
        std::vector<int> convergedPatchList_new(n_patches);

        int n_blocks = (n_patches + block_size) / block_size;

        while (activeListSize > 0) {

            // UPDATE ACTIVE PATCH
            {
                cudaMemset(converged_patch_list_dev, 0, sizeof(int) * n_patches);
                updatePatchValues3<<<activeListSize, block_size>>>(U_dev, MprimeList_dev, 7,
                                                                   activePatchList_dev, activePatchListSize_dev,
                                                                   patchElementPtr_dev, patchAdjElementIdx_dev,
                                                                   patchNodePtr_dev, patchAdjNodeIdx_dev,
                                                                   elementNodePtr_dev, elementAdjNodeIdx_dev,
                                                                   time_reduction_list_dev, time_reduction_ptr_dev,
                                                                   reduction_span, converged_patch_list_dev);

            }

            // ACTIVATE NEIGHBOURS OF CONVERGED PATCHES (MIXED)
            {
                cudaMemset(converged_patch_list_new_dev, 0, sizeof(int) * n_patches);


                activateNeighbours<<<n_blocks, block_size>>>(patchPatchPtr_dev, patchAdjPatchIdx_dev,
                                                             converged_patch_list_dev, converged_patch_list_new_dev,
                                                             n_patches);


                activeListSize = reduceActiveList(converged_patch_list_new_dev, activePatchList_dev,
                                                  activePatchListSize_dev, n_patches);
            }


            if (activeListSize == 0) break;

            //DO ONE ITERATION ON NEIGHBOURS
            {
                cudaMemset(converged_patch_list_dev, 0, sizeof(int) * n_patches);

                updatePatchValues3<<<activeListSize, block_size>>>(U_dev, MprimeList_dev, 1,
                                                                   activePatchList_dev, activePatchListSize_dev,
                                                                   patchElementPtr_dev, patchAdjElementIdx_dev,
                                                                   patchNodePtr_dev, patchAdjNodeIdx_dev,
                                                                   elementNodePtr_dev, elementAdjNodeIdx_dev,
                                                                   time_reduction_list_dev, time_reduction_ptr_dev,
                                                                   reduction_span, converged_patch_list_dev);
            }

            //RETRIEVE ACTIVE LIST (CPU)
            {
                activeListSize = reduceActiveList(converged_patch_list_dev, activePatchList_dev,
                                                  activePatchListSize_dev, n_patches);
            }

        }

        // GET FINAL VALUES
        cudaMemcpy(U.data(), U_dev, sizeof(double) * U.size(), cudaMemcpyDeviceToHost);

        // FREE ALL CUDA MEMORY
        {
            cudaFree(patchElementPtr_dev);
            cudaFree(patchAdjElementIdx_dev);
            cudaFree(patchNodePtr_dev);
            cudaFree(patchAdjNodeIdx_dev);
            cudaFree(elementNodePtr_dev);
            cudaFree(elementAdjNodeIdx_dev);
            cudaFree(patchPatchPtr_dev);
            cudaFree(patchAdjPatchIdx_dev);
            cudaFree(activePatchList_dev);
            cudaFree(U_dev);
            cudaFree(MprimeList_dev);
            cudaFree(time_reduction_ptr_dev);
            cudaFree(time_reduction_list_dev);
            cudaFree(converged_patch_list_dev);
            cudaFree(activePatchListSize_dev);
        }

        (*result) = true;
    }


    template<>
    void globalSolve<4>(const std::vector<int> &XPatches,
                        const std::vector<int> &patchElementPtr, const std::vector<int> &patchAdjElementIdx,
                        const std::vector<int> &patchNodePtr, const std::vector<int> &patchAdjNodeIdx,
                        const std::vector<int> &elementNodePtr, const std::vector<int> &elementAdjNodeIdx,
                        const std::vector<int> &patchPatchPtr, const std::vector<int> &patchAdjPatchIdx,
                        std::vector<double> &U, const std::vector<MprimeMatrix<4>> &MprimeList,
                        const std::vector<int> &time_reduction_ptr, const int reduction_span,
                        bool *result) {

        int block_size = BLOCK_SIZE;

        int n_patches = patchElementPtr.size() - 1;
        int n_elements = patchAdjElementIdx.size();
        int n_nodes = patchAdjNodeIdx.size();

        // INIT CUDA PATCH TO ELEMENT ADJACENT LIST
        int *patchElementPtr_dev, *patchAdjElementIdx_dev;
        {
            allocateAndTransfer((void **) &patchElementPtr_dev, (void *) patchElementPtr.data(),
                                sizeof(int), patchElementPtr.size());
            allocateAndTransfer((void **) &patchAdjElementIdx_dev, (void *) patchAdjElementIdx.data(),
                                sizeof(int), patchAdjElementIdx.size());
        }

        // INIT CUDA PATCH TO NODE ADJACENT LIST
        int *patchNodePtr_dev, *patchAdjNodeIdx_dev;
        {
            allocateAndTransfer((void **) &patchNodePtr_dev, (void *) patchNodePtr.data(),
                                sizeof(int), patchNodePtr.size());
            cudaMemcpy((void *) patchNodePtr.data(), patchNodePtr_dev, sizeof(int) * patchNodePtr.size(),
                       cudaMemcpyDeviceToHost);

            allocateAndTransfer((void **) &patchAdjNodeIdx_dev, (void *) patchAdjNodeIdx.data(),
                                sizeof(int), patchAdjNodeIdx.size());
        }

        // INIT CUDA ELEMENT TO NODE ADJACENT LIST
        int *elementNodePtr_dev, *elementAdjNodeIdx_dev;
        {
            allocateAndTransfer((void **) &elementNodePtr_dev, (void *) elementNodePtr.data(),
                                sizeof(int), elementNodePtr.size());
            allocateAndTransfer((void **) &elementAdjNodeIdx_dev, (void *) elementAdjNodeIdx.data(),
                                sizeof(int), elementAdjNodeIdx.size());
        }

        // INIT CUDA PATCH TO PATCH ADJACENT LIST
        int *patchPatchPtr_dev, *patchAdjPatchIdx_dev;
        {
            allocateAndTransfer((void **) &patchPatchPtr_dev, (void *) patchPatchPtr.data(),
                                sizeof(int), patchPatchPtr.size());
            allocateAndTransfer((void **) &patchAdjPatchIdx_dev, (void *) patchAdjPatchIdx.data(),
                                sizeof(int), patchAdjPatchIdx.size());
        }

        // INIT CUDA MPRIMEMATRIX PER ELEMENT LIST
        MprimeMatrix<4> *MprimeList_dev;
        {
            allocateAndTransfer((void **) &MprimeList_dev, (void *) MprimeList.data(),
                                sizeof(MprimeMatrix<4>), MprimeList.size());
        }

        // INIT CUDA TIME REDUCTION LIST AND POINTER LIST
        int *time_reduction_ptr_dev;
        double *time_reduction_list_dev;
        {
            allocateAndTransfer((void **) &time_reduction_ptr_dev, (void *) time_reduction_ptr.data(),
                                sizeof(int), time_reduction_ptr.size());

            // init all time reduction list at MAXF
            int list_size = reduction_span * n_nodes;
            cudaMalloc((void **) &time_reduction_list_dev, sizeof(double) * list_size);
            int n_blocks = (list_size + block_size) / block_size;
            initTimeReductionList<<<n_blocks, block_size>>>(time_reduction_list_dev, list_size);
        }

        // INIT CUDA SUPPORT LISTS FOR CONVERGED NODES AND PATCHES
        int *converged_patch_list_dev, *converged_patch_list_new_dev;
        {
            cudaMalloc((void **) &converged_patch_list_dev, sizeof(int) * n_patches);
            cudaMalloc((void **) &converged_patch_list_new_dev, sizeof(int) * n_patches);
        }

        // INIT CUDA TIME
        double *U_dev;
        {
            allocateAndTransfer((void **) &U_dev, (void *) U.data(),
                                sizeof(double), U.size());
        }

        // INIT CUDA ACTIVE LIST AND SIZE
        std::vector<int> activePatchList(n_patches);
        std::copy(XPatches.begin(), XPatches.end(), activePatchList.data());

        int activeListSize = XPatches.size();
        int *activePatchList_dev, *activePatchListSize_dev;
        {
            allocateAndTransfer((void **) &activePatchList_dev, (void *) activePatchList.data(),
                                sizeof(int), n_patches);
            allocateAndTransfer((void **) &activePatchListSize_dev, (void *) &activeListSize,
                                sizeof(int), 1);
        }


        std::vector<int> convergedPatchList(n_patches);
        std::vector<int> convergedPatchList_new(n_patches);

        int n_blocks = (n_patches + block_size) / block_size;

        while (activeListSize > 0) {

            // UPDATE ACTIVE PATCH
            {
                cudaMemset(converged_patch_list_dev, 0, sizeof(int) * n_patches);
                updatePatchValues4<<<activeListSize, block_size>>>(U_dev, MprimeList_dev, 7,
                                                                   activePatchList_dev, activePatchListSize_dev,
                                                                   patchElementPtr_dev, patchAdjElementIdx_dev,
                                                                   patchNodePtr_dev, patchAdjNodeIdx_dev,
                                                                   elementNodePtr_dev, elementAdjNodeIdx_dev,
                                                                   time_reduction_list_dev, time_reduction_ptr_dev,
                                                                   reduction_span, converged_patch_list_dev);

            }

            // ACTIVATE NEIGHBOURS OF CONVERGED PATCHES (GPU)
            {
                cudaMemset(converged_patch_list_new_dev, 0, sizeof(int) * n_patches);
                activateNeighbours<<<n_blocks, block_size>>>(patchPatchPtr_dev, patchAdjPatchIdx_dev,
                                                             converged_patch_list_dev, converged_patch_list_new_dev,
                                                             n_patches);

                activeListSize = reduceActiveList(converged_patch_list_new_dev, activePatchList_dev,
                                                  activePatchListSize_dev, n_patches);

            }


            if (activeListSize == 0) break;

            //DO ONE ITERATION ON NEIGHBOURS
            {
                cudaMemset(converged_patch_list_dev, 0, sizeof(int) * n_patches);

                updatePatchValues4<<<activeListSize, block_size>>>(U_dev, MprimeList_dev, 1,
                                                                   activePatchList_dev, activePatchListSize_dev,
                                                                   patchElementPtr_dev, patchAdjElementIdx_dev,
                                                                   patchNodePtr_dev, patchAdjNodeIdx_dev,
                                                                   elementNodePtr_dev, elementAdjNodeIdx_dev,
                                                                   time_reduction_list_dev, time_reduction_ptr_dev,
                                                                   reduction_span, converged_patch_list_dev);
            }

            //RETRIEVE ACTIVE LIST (gpU)
            {
                activeListSize = reduceActiveList(converged_patch_list_dev, activePatchList_dev,
                                                  activePatchListSize_dev, n_patches);
            }

        }

        // GET FINAL VALUES
        cudaMemcpy(U.data(), U_dev, sizeof(double) * U.size(), cudaMemcpyDeviceToHost);

        // FREE ALL CUDA MEMORY
        {
            cudaFree(patchElementPtr_dev);
            cudaFree(patchAdjElementIdx_dev);
            cudaFree(patchNodePtr_dev);
            cudaFree(patchAdjNodeIdx_dev);
            cudaFree(elementNodePtr_dev);
            cudaFree(elementAdjNodeIdx_dev);
            cudaFree(patchPatchPtr_dev);
            cudaFree(patchAdjPatchIdx_dev);
            cudaFree(activePatchList_dev);
            cudaFree(U_dev);
            cudaFree(MprimeList_dev);
            cudaFree(time_reduction_ptr_dev);
            cudaFree(time_reduction_list_dev);
            cudaFree(converged_patch_list_dev);
            cudaFree(activePatchListSize_dev);
        }

        (*result) = true;
    }


}
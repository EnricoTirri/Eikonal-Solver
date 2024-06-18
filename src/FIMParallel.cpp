
#if FIMP

#include "EikonalSolver.hpp"
#include <iostream>
#include "EikonalTraits.hpp"
#include "EikonalHeapComparator.cpp"
#include "LocalSolver.hpp"
#include <cmath>
#include <chrono>
namespace Eikonal {
    using Point = Eikonal::Traits::Point;

    template<int MESH_SIZE>
    void EikonalSolver<MESH_SIZE>::print_spec() {
        std::cout << "Fast Iterative Method OpenMP parallelized Eikonal solver" << std::endl;
    }

    template<int MESH_SIZE>
    inline void
    initialize(std::vector<double> &U, const std::vector<int> &X, const Mesh<MESH_SIZE> &data, std::vector<int> &active, const double &maxf) {
            // Init values
#pragma omp for
            for (int i = 0; i < U.size(); ++i)
                U[i] = maxf;

            // Init actives
#pragma omp for
            for (int i = 0; i < active.size(); ++i)
                active[i] = 0;

            // Activate starting agglomerates
#pragma omp for
            for (int i = 0; i < X.size(); ++i) {
                int pId = X[i];
                U[pId] = 0;
                int elRS = data.adjPointPtr[pId];
                int elRE = data.adjPointPtr[pId + 1];
//#pragma omp parallel for
                for (int j = elRS; j < elRE; ++j) {
                    int eId = data.pointAdjacentElementList[j];
                    int pRS = data.adjElementPtr[eId];
                    int pRE = data.adjElementPtr[eId + 1];
                    for (int k = pRS; k < pRE; ++k) {
                        int tpId = data.elementAdjacentPointList[k];
                        active[tpId] = 1;
                    }
                }
            }
    }

    inline void applyScanAndPack(std::vector<int> &active, std::vector<int> &reducedActive) {
        //PARALLEL SCAN FOR ACTIVE LIST PACK INDICES
            for (int i = 0; i < std::log2(active.size()); ++i) {
                int stride = static_cast<int>(pow(2, i));
                int step = 2 * stride;
#pragma omp for simd
                for (int k = step - 1; k < active.size(); k += step) {
                    active[k] += active[k - stride];
                }
            }

            for (int i = std::log2(active.size()); i > 0; i--) {
                int stride = static_cast<int>(pow(2, i));
                int step = stride / 2;
#pragma omp for simd
                for (int k = stride - 1; k < active.size() - step; k += stride) {
                    active[k + step] += active[k];
                }
            }
#pragma omp single
            {
                reducedActive.resize(active.back());
            }
#pragma omp barrier

#pragma omp for
            for (int i = 0;i < active.size();++i) {
                if (active[i] > 0) {
                    if (i == 0 || active[i - 1] < active[i]) {
                        reducedActive[active[i] - 1] = i;
                    }
                }
            }
    }

    template<int MESH_SIZE>
    inline void
    updateAgglomerate(const int &agglID, std::vector<double> &newValues, const std::vector<double> &U, const Mesh<MESH_SIZE> &data, const int &iter) {
        int v = agglID;

        int elRangeStart = data.adjPointPtr[v];
        int elRangeEnd = data.adjPointPtr[v + 1];

        newValues.resize(elRangeEnd - elRangeStart);


            for (int i = elRangeStart; i < elRangeEnd; i++) {

                int elID = data.pointAdjacentElementList[i];

                int pntRangeStart = data.adjElementPtr[elID];

                std::array<Point, MESH_SIZE> base;
                int k = 0;
                std::array<double, MESH_SIZE> values;
                for (int j = 0; j < MESH_SIZE; j++) {
                    int tpId = data.elementAdjacentPointList[pntRangeStart + j];
                    if (tpId != v) {
                        base[k] = data.points[tpId];
                        values[k] = U[tpId];
                        k++;
                    }
                }
                base[MESH_SIZE - 1] = data.points[v];

                Traits::VelocityM M;
                M << 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;

                LocalSolver<MESH_SIZE> solver(M, base, values, iter, 10e-9);
                double sol = solver();
                newValues[i - elRangeStart] = sol;
            }

            for (int step = 1; step < newValues.size(); step *= 2) {
                for (int j = 0; j < newValues.size() - step; j += step) {
                    if (newValues[j] > newValues[j + step])
                        newValues[j] = newValues[j + step];
                }
            }

    }


    template<int MESH_SIZE>
    bool EikonalSolver<MESH_SIZE>::solve(std::vector<double> &U, const std::vector<int> &X,
                                         const Mesh<MESH_SIZE> &data) {


        auto start = std::chrono::high_resolution_clock::now();
        // Check if all points belong to mesh
        for (auto &point: X) {
            if (point >= data.points.size()) {
                printf("error on initial point id: %d does not belong to mesh\n", point);
                return false;
            }
        }

        std::vector<int> active;
        std::vector<int> reducedActive;
        std::vector<double> converged;

        // ---------------------- INITIALIZATION --------------------- //
        U.resize(data.points.size());
        active.resize(data.points.size());
        converged.resize(data.points.size());

        auto end = std::chrono::high_resolution_clock::now();
        this->prepare = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        start = std::chrono::high_resolution_clock::now();
#pragma omp parallel
        {
        initialize<MESH_SIZE>(U, X, data, active, MAXF);
        // ----------------------------------------------------------- //

        // ------------------- PACK OF ACTIVE LIST ------------------- //
            applyScanAndPack(active, reducedActive);
        // ----------------------------------------------------------- //


#ifdef SOLVER_VERBOSE
        int iteration = 0;
#endif

        while (!reducedActive.empty()) {

#ifdef SOLVER_VERBOSE
            iteration++;
            std::cout << "Iteration: " << iteration << "\tActiveList size: " << active.size() << std::endl;
#endif

                // --------------- UPDATE ALL ACTIVE AGGLOMERATES ----------------- //
#pragma omp for simd
                for (int i = 0; i < active.size(); ++i) {
                    converged[i] = -1;
                    active[i] = 0;
                }

#pragma omp for
                for (int i = 0; i < reducedActive.size(); ++i) {
                    std::vector<double> newValues;

                    int aggID = reducedActive[i];

                    updateAgglomerate<MESH_SIZE>(aggID, newValues, U, data, 5000);

                    if (U[aggID] > newValues[0]) {
                        converged[aggID] = newValues[0];
                    }
                }
                // ---------------------------------------------------------------- //

                // - IF AGGLOMERATE HAS CONVERGED, ACTIVATE ADJACENT AGGLOMERATES - //

#pragma omp for
                for (int i = 0; i < reducedActive.size(); ++i) {
                    int pId = reducedActive[i];
                    if (converged[pId] > 0) {
                        U[pId] = converged[pId];

                        int elRS = data.adjPointPtr[pId];
                        int elRE = data.adjPointPtr[pId + 1];
                        for (int j = elRS; j < elRE; ++j) {
                            int eId = data.pointAdjacentElementList[j];
                            int pRS = data.adjElementPtr[eId];
                            int pRE = data.adjElementPtr[eId + 1];
                            for (int k = pRS; k < pRE; ++k) {
                                int tpId = data.elementAdjacentPointList[k];
                                active[tpId] = 1;
                            }
                        }
                    }
                }

#pragma omp for simd
                for (int i = 0; i < reducedActive.size(); ++i) {
                    int pId = reducedActive[i];
                    if(converged[i] > 0)
                        active[pId] = 0;
                }

            // ------------------- PACK OF ACTIVE LIST ------------------- //
            applyScanAndPack(active, reducedActive);
            // ----------------------------------------------------------- //


            // -------- CHECK IF ACTIVE AGGLOMERATES CONVERGE ----------------- //
#pragma omp for simd
                for (int i = 0; i < active.size(); i++) {
                    converged[i] = -1;
                    active[i] = 0;
                }

#pragma omp for
                for (int i = 0; i < reducedActive.size(); ++i) {
                    std::vector<double> newValues;

                    int aggID = reducedActive[i];

                    updateAgglomerate<MESH_SIZE>(aggID, newValues, U, data, 1);

                    if (U[aggID] > newValues[0]) {
                        converged[aggID] = newValues[0];
                    }
                }

#pragma omp for simd
                for (int i = 0; i < reducedActive.size(); ++i) {
                    int ptID = reducedActive[i];
                    if (converged[ptID] > 0)
                        active[ptID] = 1;
                }
            // ---------------------------------------------------------------- //

            // ----------------------------------------------------------------- //

            // ------------------- PACK OF ACTIVE LIST ------------------- //
            applyScanAndPack(active, reducedActive);
            // ----------------------------------------------------------- //

//#pragma omp barrier

            }
        }
        end = std::chrono::high_resolution_clock::now();
        this->compute = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        return true;
    }

}
#endif


#if FMMO

#include "EikonalSolver.hpp"
#include <queue>
#include <iostream>
#include "EikonalTraits.hpp"
#include "EikonalHeapComparator.cpp"
#include "OptimizedLocalSolver.hpp"

#include <chrono>
namespace Eikonal {
    template<int MESH_SIZE>
    void EikonalSolver<MESH_SIZE>::print_spec() {
        std::cout << "Fast Marching Method Eikonal solver, Local Optimized" << std::endl;
    }

    template<int MESH_SIZE>
    bool EikonalSolver<MESH_SIZE>::solve(std::vector<double> &U, const std::vector<int> &X,
                                         const Mesh<MESH_SIZE> &data) {


        using Point = typename Eikonal::Traits::Point;


        //init start prepare time
        auto start = std::chrono::high_resolution_clock::now();

        // Check if pointId belongs to mesh
        for (auto &pointId: X) {
            if (pointId >= data.points.size()) {
                printf("error on initial point id: %d does not belong to mesh\n", pointId);
                return false;
            }
        }

        // Initialize minHeap
        U.resize(data.points.size());
        std::priority_queue<int, std::vector<int>, EikonalHeapComparator> minHeap((EikonalHeapComparator(U)));
        std::fill(U.begin(), U.end(), MAXF);

        std::vector<OptimizedLocalSolver<MESH_SIZE>> solvers;

        Traits::VelocityM M;
        M << 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;

        std::array<Traits::Point, MESH_SIZE> points;
        for (int i = 0; i < data.adjElementPtr.size() - 1; ++i) {
            int ptrStart = data.adjElementPtr[i];
            int ptrEnd = data.adjElementPtr[i + 1];

            for (int j = ptrStart; j < ptrEnd; ++j) {
                points[j - ptrStart] = data.points[data.elementAdjacentPointList[j]];
            }

            solvers.push_back(OptimizedLocalSolver<MESH_SIZE>{5000, 10e-6, M, points});
        }


        for (const auto &i: X) {
            U[i] = 0;
            minHeap.push(i);
        }
        auto end = std::chrono::high_resolution_clock::now();

        this->prepare = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        start = std::chrono::high_resolution_clock::now();
#ifdef SOLVER_VERBOSE
        int iteration = 0;
#endif

        std::vector<bool> L_set(data.points.size(), false);
        std::vector<bool> L_in(data.points.size(), false);

        while (!minHeap.empty()) {

#ifdef SOLVER_VERBOSE
            iteration++;
            std::cout << "Iteration: " << iteration << "\tActiveList size: " << minHeap.size() << std::endl;
#endif


            // Take the index of the not already explored point with lower time
            int minPointID = minHeap.top();
            minHeap.pop();
            L_set[minPointID] = true;


            // Add all elements the point belongs to, to the elements that must be updated
            std::vector<int> neighbors;
            int elRangeStart = data.adjPointPtr[minPointID];
            int elRangeEnd = data.adjPointPtr[minPointID + 1];

            for (std::size_t j = elRangeStart; j < elRangeEnd; j++) {
                neighbors.push_back(data.pointAdjacentElementList[j]);
            }

            for (const auto &eleId: neighbors) {
                // Take range of pntId of the element
                int pointRangeStart = data.adjElementPtr[eleId];
                int pointRangeEnd = data.adjElementPtr[eleId + 1];

                OptimizedLocalSolver<MESH_SIZE> solver = solvers[eleId];

                // Otherwise solve local problem with this point as unknown and the others as base
                std::array<double, MESH_SIZE> values;
                for (std::size_t j = 0; j < MESH_SIZE; j++) {
                    int tpId = data.elementAdjacentPointList[pointRangeStart + j];
                    values[j] = U[tpId];
                }

                for (int pi = pointRangeStart; pi < pointRangeEnd; ++pi) {
                    int pointID = data.elementAdjacentPointList[pi];

                    // If point of the element has been already explored then continue
                    if (L_set[pointID]) continue;

                    int k = pi - pointRangeStart;

                    double sol = solver(k, values);

                    // If solution does not improve then skip
                    if (sol > U[pointID]) continue;

                    // Otherwise update solution
                    U[pointID] = sol;
                    values[k] = sol;

                    //
                    if (!L_in[pointID]) {
                        L_in[pointID] = true;
                        minHeap.push(pointID);
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
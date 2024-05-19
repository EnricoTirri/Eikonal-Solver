
#if FMM

#include "EikonalSolver.hpp"
#include <queue>
#include <iostream>
#include "EikonalTraits.hpp"
#include "EikonalHeapComparator.cpp"
#include "LocalSolver.hpp"

namespace Eikonal {
    template<int MESH_SIZE>
    void EikonalSolver<MESH_SIZE>::print_spec() {
        std::cout << "Fast Marching Method Eikonal solver" << std::endl;
    }

    template<int MESH_SIZE>
    bool EikonalSolver<MESH_SIZE>::solve(std::vector<double> &U, const std::vector<int> &X,
                                              const Mesh<MESH_SIZE> &data) {

        using Point = typename Eikonal::Traits::Point;

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
        for (const auto &i: X) {
            U[i] = 0;
            minHeap.push(i);
        }

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
            int elRangeEnd = data.adjPointPtr[minPointID+1];

            for (std::size_t j = elRangeStart; j < elRangeEnd; j++) {
                neighbors.push_back(data.pointAdjacentElementList[j]);
            }

            for (const auto &eleId: neighbors) {
                // Take range of pntId of the element
                int pointRangeStart = data.adjElementPtr[eleId];
                int pointRangeEnd = data.adjElementPtr[eleId + 1];

                for (int pi = pointRangeStart; pi < pointRangeEnd; ++pi) {
                    int pointID = data.elementAdjacentPointList[pi];

                    // If point of the element has been already explored then continue
                    if (L_set[pointID]) continue;

                    // Otherwise solve local problem with this point as unknown and the others as base
                    std::array<Point, MESH_SIZE> base;
                    std::size_t k = 0;
                    Eigen::Matrix<double, MESH_SIZE, 1> values;
                    for (std::size_t j = 0; j < MESH_SIZE; j++) {
                        int tpId = data.elementAdjacentPointList[pointRangeStart+j];
                        if (tpId == pointID) {
                            base[MESH_SIZE - 1] = data.points[pointID];
                        } else {
                            base[k] = data.points[tpId];
                            values[k] = U[tpId];
                            k++;
                        }
                    }

                    Traits::VelocityM M;
                    M << 1.0, 0.0, 0.0,
                            0.0, 1.0, 0.0,
                            0.0, 0.0, 1.0;

                    LocalSolver<MESH_SIZE> solver(M, base, values, 5000, 10e-6);

                    double sol = solver();

                    // If solution does not improve then skip
                    if (sol > U[pointID]) continue;

                    // Otherwise update solution
                    U[pointID] = sol;

                    //
                    if (!L_in[pointID]) {
                        L_in[pointID] = true;
                        minHeap.push(pointID);
                    }
                }
            }
        }
        return true;
    }
}
#endif
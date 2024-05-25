
#if FMMP

#include "EikonalSolver.hpp"
#include <queue>
#include <iostream>
#include "EikonalTraits.hpp"
#include "EikonalHeapComparator.cpp"
#include "LocalSolver.hpp"

namespace Eikonal {

    class PointElements {
    public:
        size_t punto;
        std::vector<int> elements_legacy;

        PointElements(size_t p) : punto(p) {}
    };

    template<int MESH_SIZE>
    void EikonalSolver<MESH_SIZE>::print_spec() {
        std::cout << "Fast Marching Omp-Parallel Method Eikonal solver" << std::endl;
    }

    template<int MESH_SIZE>
    bool EikonalSolver<MESH_SIZE>::solve(std::vector<double> &U,
                                         const std::vector<int> &X,
                                         const Mesh<MESH_SIZE> &data) {

        typedef typename Traits::Point Point;

        // Check if pointId belongs to mesh
        for (auto &point: X) {
            if (point >= data.points.size()) {
                printf("error on initial point id: %d does not belong to mesh\n", point);
                return false;
            }
        }

        // Initialize minHeap
        U.reserve(data.points.size());
        std::priority_queue<int, std::vector<int>, EikonalHeapComparator> minHeap((EikonalHeapComparator(U)));
        std::fill(U.begin(), U.end(), MAXF);
        for (const auto &i: X) {
            U[i] = 0;
            minHeap.push(i);
        }

        // Initialize progress info
        int step_count = data.adjPointPtr.size() / 50;
        int count = 0;
        int step = 0;

        std::vector<bool> L_set(data.points.size(), false);
        std::vector<bool> L_in(data.points.size(), false);
        while (!minHeap.empty()) {
            int minPointID = minHeap.top();
            minHeap.pop();
            L_set[minPointID] = true;

            // Update progress
            ++count;
            if (step_count == count - 1) {
                count = 0;
                step += 2;
                std::cout << '\r';
                std::cout << step << "% ";
                std::cout.flush();
            }

            std::vector<int> neighbors;
            int elRangeStart = data.adjPointPtr[minPointID];
            int elRangeEnd = data.adjPointPtr[minPointID+1];

            for (std::size_t j = elRangeStart; j < elRangeEnd; j++) {
                neighbors.push_back(data.pointAdjacentElementList[j]);
            }

            std::vector<PointElements> activepoints;

            // Find distinct point that have to be updated
            for (size_t index1 = 0; index1 < neighbors.size(); ++index1) {
                int elID = neighbors[index1];
                int pntRangeStart = data.adjElementPtr[elID];

                for (size_t index2 = 0; index2 < MESH_SIZE; ++index2) {
                    int pointID = data.elementAdjacentPointList[pntRangeStart + index2];

                    if (!L_in[pointID]) {
                        if (activepoints.empty()) {
                            PointElements activepoint(pointID);
                            activepoint.elements_legacy.push_back(elID);
                            activepoints.push_back(activepoint);
                        }

                        int inserted = 0;
                        for (auto & activepoint : activepoints) {
                            if (activepoint.punto == pointID) {
                                activepoint.elements_legacy.push_back(elID);
                                ++inserted;
                            }
                        }

                        if (!inserted && !L_set[pointID]) {
                            PointElements activepoint(pointID);
                            activepoint.elements_legacy.push_back(elID);
                            activepoints.push_back(activepoint);
                        }
                    }
                }
            }


#pragma omp parallel for schedule(static) num_threads(N_THREADS)
            for (auto & activepoint : activepoints) {
                const auto pointID = activepoint.punto;
                for (int elID : activepoint.elements_legacy) {
                    int pointRangeStart = data.adjElementPtr[elID];

                    std::array<Point, MESH_SIZE> base;
                    std::size_t k = 0;
                    std::array<double, MESH_SIZE> values;
                    for (std::size_t j = 0; j < MESH_SIZE; j++) {
                        int tpId = data.elementAdjacentPointList[pointRangeStart+j];
                        if (tpId != pointID) {
                            base[k] = data.points[tpId];
                            values[k] = U[tpId];
                            k++;
                        }
                    }
                    base[MESH_SIZE - 1] = data.points[pointID];

                    Traits::VelocityM M;
                    M << 1.0, 0.0, 0.0,
                            0.0, 1.0, 0.0,
                            0.0, 0.0, 1.0;

                    LocalSolver <MESH_SIZE> solver(M, base, values, 5000, 10e-6);

                    auto sol = solver();

                    if (sol > U[pointID]) continue;
                    U[pointID] = sol;

                }
            }

            for (auto &activepoint: activepoints) {
                if (!L_in[activepoint.punto]) {
                    L_in[activepoint.punto] = true;
                    minHeap.push(activepoint.punto);
                }
            }

        }

        std::cout << std::endl;
        return true;
    }
}
#endif

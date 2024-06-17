
#if FIM

#include "EikonalSolver.hpp"
#include <iostream>
#include "EikonalTraits.hpp"
#include "EikonalHeapComparator.cpp"
#include "LocalSolver.hpp"

#include <chrono>
namespace Eikonal {
    using Point = typename Eikonal::Traits::Point;

    template<int MESH_SIZE>
    void EikonalSolver<MESH_SIZE>::print_spec() {
        std::cout << "Fast Iterative Method Eikonal solver" << std::endl;
    }

    void removeDuplicates(std::vector<int> &v) {
        sort(v.begin(), v.end());
        v.resize(distance(v.begin(), unique(v.begin(), v.end())));
    }

    template<int MESH_SIZE>
    bool EikonalSolver<MESH_SIZE>::solve(std::vector<double> &U, const std::vector<int> &X,
                                         const Mesh<MESH_SIZE> &data) {



        auto start = std::chrono::high_resolution_clock::now();
        for (auto &point: X) {
            if (point >= data.points.size()) {
                printf("error on initial point id: %d does not belong to mesh\n", point);
                return false;
            }
        }

        U.resize(data.points.size());
        std::fill(U.begin(), U.end(), MAXF);

        auto end = std::chrono::high_resolution_clock::now();
       this->prepare = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        start = std::chrono::high_resolution_clock::now();

        std::vector<int> active;


        for (int pntID: X) {
            U[pntID] = 0;

            int elRangeStart = data.adjPointPtr[pntID];
            int elRangeEnd = data.adjPointPtr[pntID + 1];

            // for each neighbour
            for (std::size_t j = elRangeStart; j < elRangeEnd; j++) {
                int elID = data.pointAdjacentElementList[j];

                // find points of element
                int pntRangeStart = data.adjElementPtr[elID];
                int pntRangeEnd = data.adjElementPtr[elID + 1];

                // for each point
                for (int i = pntRangeStart; i < pntRangeEnd; ++i) {
                    int tpntID = data.elementAdjacentPointList[i];

                    // add to active set
                    if(tpntID!=pntID)
                        active.push_back(tpntID);
                }
            }
        }

        removeDuplicates(active);

        std::vector<int> activeNew;
        std::vector<int> converged;

#ifdef SOLVER_VERBOSE
        int iteration = 0;
#endif

        std::vector<bool> activeSet(data.points.size());
        while (!active.empty()) {
            std::fill(activeSet.begin(), activeSet.end(), false);


#ifdef SOLVER_VERBOSE
            iteration++;
            std::cout << "Iteration: " << iteration << "\tActiveList size: " << active.size() << std::endl;
#endif

            converged.clear();
            for (const int &v: active) {
                double min = U[v];

                int elRangeStart = data.adjPointPtr[v];
                int elRangeEnd = data.adjPointPtr[v + 1];
                for (std::size_t i = elRangeStart; i < elRangeEnd; i++) {
                    int elID = data.pointAdjacentElementList[i];

                    int pntRangeStart = data.adjElementPtr[elID];

                    std::array<Point, MESH_SIZE> base;
                    std::size_t k = 0;
                    std::array<double, MESH_SIZE> values;
                    for (std::size_t j = 0; j < MESH_SIZE; j++) {
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

                    LocalSolver<MESH_SIZE> solver(M, base, values, 5000, 10e-6);
                    double sol = solver();
                    if (sol > min) continue;
                    min = sol;
                }

                if (min < U[v]) {
                    U[v] = min;
                    converged.push_back(v);
                }
            }

            activeNew.clear();
            for (const int &v: converged) {
                int elRangeStart = data.adjPointPtr[v];
                int elRangeEnd = data.adjPointPtr[v + 1];

                // for each neighbour
                for (std::size_t j = elRangeStart; j < elRangeEnd; j++) {
                    int elID = data.pointAdjacentElementList[j];

                    // find points of element
                    int pntRangeStart = data.adjElementPtr[elID];
                    int pntRangeEnd = data.adjElementPtr[elID + 1];

                    // for each point
                    for (int i = pntRangeStart; i < pntRangeEnd; ++i) {
                        int tpntID = data.elementAdjacentPointList[i];

                        if (!activeSet[tpntID]) {
                            activeSet[tpntID] = true;
                            // add to active set
                            activeNew.push_back(tpntID);
                        }
                    }
                }
            }
            //removeDuplicates(activeNew);

            active.clear();
            for (const int &v: activeNew) {
                double min = U[v];

                int elRangeStart = data.adjPointPtr[v];
                int elRangeEnd = data.adjPointPtr[v + 1];

                for (std::size_t i = elRangeStart; i < elRangeEnd; i++) {
                    int elID = data.pointAdjacentElementList[i];

                    int pntRangeStart = data.adjElementPtr[elID];

                    std::array<Point, MESH_SIZE> base;
                    std::size_t k = 0;
                    std::array<double, MESH_SIZE> values;
                    for (std::size_t j = 0; j < MESH_SIZE; j++) {
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

                    LocalSolver<MESH_SIZE> solver(M, base, values, 1, 10e-6);
                    double sol = solver();

                    if (sol < min)
                        min = sol;
                }

                if (min < U[v]) {
                    active.push_back(v);
                }
            }
        }

            end = std::chrono::high_resolution_clock::now();
       this->compute = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        return true;
    }

}

#endif

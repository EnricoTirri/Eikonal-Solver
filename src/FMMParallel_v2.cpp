//
// Created by Enrico on 19/12/2023.
//

#if FMMP2

#include "EikonalSolver.hpp"
#include <queue>
#include <iostream>
#include "Eikonal_traits.hpp"
#include "SimplexData.hpp"
#include "solveEikonalLocalProblem.hpp"
#include <execution>
#include "EikonalHeapComparator.cpp"

class ParallelStruct {
public:
    size_t point;
    double value;

    // Constructor with default values
    ParallelStruct(size_t p = 0, double v = 0.0) : point(p), value(v) {}
};


class PointElements {
public:
    size_t punto;
    std::vector<int> elements;

    PointElements(size_t p) : punto(p) {}
};

template<int DIM, int MESH_SIZE>
void EikonalSolver<DIM, MESH_SIZE>::print_spec(){
    std::cout << "Fast Marching Omp-Parallel Method (v2) Eikonal solver" << std::endl;
}

template<int DIM, int MESH_SIZE>
bool EikonalSolver<DIM,MESH_SIZE>::solve(std::vector<double> &U,
                         const std::vector<int>& X,
                         const Mesh<DIM, MESH_SIZE> &data) {

    typedef typename Eikonal_traits<DIM>::Point Point;
    for (auto &point: X) {
        if (point >= data.index.size()) {
            printf("error on initial point: %f %f does not belong to mesh\n", data.points[point].x(),
                   data.points[point].y());
            return false;
        }
    }

    U.reserve(data.index.size());
    std::priority_queue<int, std::vector<int>, EikonalHeapComparator>
            minHeap((EikonalHeapComparator(U)));
    std::fill(std::execution::seq, U.begin(), U.end(), MAXF);
    for (const auto &i: X) {
        U[i] = 0;
        minHeap.push(i);
    }

    int step_count = data.index.size() / 50;
    int count = 0;
    int step = 0;

    std::vector<bool> L_set(data.points.size(), false);
    std::vector<bool> L_in(data.points.size(), false);
    while (!minHeap.empty()) {
        int i = minHeap.top();
        minHeap.pop();
        L_set[i] = true;
        ++count;
        //   std::cout << count << " " << i << std::endl;
        if (step_count == count - 1) {
            count = 0;
            step += 2;
            std::cout << '\r';
            std::cout << step << "% ";
            std::cout.flush();
        }

        std::vector<int> neighbors;
        std::size_t start, end;
        start = data.index.at(i).start;
        end = data.index.at(i).end;
        for (std::size_t j = start; j < end; j++) {
            neighbors.push_back(data.adjacentList[j]);
        }
        bool error = false;

        std::vector<PointElements> activepoints;

        for (size_t index1 = 0; index1 < neighbors.size(); ++index1) {
            for (size_t index2 = 0; index2 < MESH_SIZE; ++index2) {
                if (!L_in[data.elements[neighbors[index1]][index2]]) {
                    if (activepoints.empty()) {
                        PointElements activepoint(data.elements[neighbors[index1]][index2]);
                        activepoint.elements.push_back(neighbors[index1]);
                        activepoints.push_back(activepoint);
                    }
                    int inserted = 0;
                    for (size_t index3 = 0; index3 < activepoints.size(); ++index3) {
                        if (activepoints[index3].punto == data.elements[neighbors[index1]][index2]) {
                            activepoints[index3].elements.push_back(neighbors[index1]);
                            ++inserted;
                        }
                    }
                    if (!inserted && !L_set[data.elements[neighbors[index1]][index2]]) {
                        PointElements activepoint(data.elements[neighbors[index1]][index2]);
                        activepoint.elements.push_back(neighbors[index1]);
                        activepoints.push_back(activepoint);
                    }
                }
            }
        }

#pragma omp parallel for schedule(static) num_threads(N_THREADS)
        for (size_t index1 = 0; index1 < activepoints.size(); ++index1) {
            const auto point = activepoints[index1].punto;
            for (size_t index2 = 0; index2 < activepoints[index1].elements.size(); ++index2) {
                const auto m_element = activepoints[index1].elements[index2];

                // if point in L continue
                // solve local problem with this point as unknown and the others as base
                std::array<Point, MESH_SIZE> base;
                std::size_t k = 0;
                Eigen::Matrix<double, MESH_SIZE, 1> values;
                for (std::size_t j = 0; j < MESH_SIZE; j++) {
                    if ((data.elements[m_element])[j] != point) {
                        base[k] = data.points[(data.elements[m_element])[j]];
                        values[k] = U[(data.elements[m_element])[j]];
                        k++;
                    }
                }
                base[MESH_SIZE - 1] = data.points[point];
                typename Eikonal::Eikonal_traits<DIM>::MMatrix M;
                if constexpr (DIM == 2)
                    M << 1.0, 0.0,
                            0.0, 1.;
                else if constexpr (DIM == 3)
                    M << 1.0, 0.0, 0.0,
                            0.0, 1.0, 0.0,
                            0.0, 0.0, 1.0;
                Eikonal::SimplexData<DIM, MESH_SIZE> simplex{base, M};
                Eikonal::solveEikonalLocalProblem<DIM, MESH_SIZE> solver{std::move(simplex),
                                                                         values};

                auto sol = solver();
                // if no descent direction or no convergence kill the process
                if (sol.status != 0) {
                    printf("error on convergence\n");
                    error = true;
                    // continue;
                }
                if (sol.value > U[point]) continue;
                U[point] = sol.value;

            }
        }
        if (error)
            return false;

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

#endif

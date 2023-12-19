//
// Created by Enrico on 19/12/2023.
//

#if FMMP

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
    std::cout << "Fast Marching Omp-Parallel Method Eikonal solver" << std::endl;
}

template<int DIM, int MESH_SIZE>
bool EikonalSolver<DIM,MESH_SIZE>::solve(std::vector<double> &U,
                        const std::vector<int>& X,
                        const Mesh<DIM, MESH_SIZE> &data) {

    using Point = typename Eikonal_traits<DIM>::Point;
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
    for (int i = 0; i < data.points.size(); i++) {
        U[i] = MAXF;
    }
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
        if (step_count == count - 1) {
            count = 0;
            step += 2;
            std::cout << '\r';
            std::cout << step << "% ";
            std::cout.flush();
        }

        // take time value of point p
        // find neighbors of L[i] and get the base (the DIMENSION points with the smallest value of U
        std::vector<int> neighbors;
        std::size_t start, end;
        start = data.index.at(i).start;
        end = data.index.at(i).end;
        for (std::size_t j = start; j < end; j++) {
            neighbors.push_back(data.adjacentList[j]);
        }
        bool error = false;
        // maybe we can parallelize this

        ParallelStruct parallelarray[neighbors.size()][MESH_SIZE];

#pragma omp parallel for collapse(2) num_threads(8)
        for (size_t index1 = 0; index1 < neighbors.size(); ++index1) {
            for (size_t index2 = 0; index2 < MESH_SIZE; ++index2) {
                const auto point = data.elements[neighbors[index1]][index2];
                const auto m_element = neighbors[index1];
                ParallelStruct p(point, MAXF);
                parallelarray[index1][index2] = p;

                if (point == i || L_set[point])
                    continue;
                // if point in L continue
                // solve local problem with this point as unknown and the others as base
                std::array<Point, MESH_SIZE> base;
                std::size_t k = 0;
                Eigen::Matrix<double, MESH_SIZE, 1> values;
                for (std::size_t j = 0; j < MESH_SIZE; j++) {
                    if ((data.elements[m_element])[j] == point) {
                        base[MESH_SIZE - 1] = data.points[point];
                    } else {
                        base[k] = data.points[(data.elements[m_element])[j]];
                        values[k] = U[(data.elements[m_element])[j]];
                        k++;
                    }
                }
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
                p.value = sol.value;

                parallelarray[index1][index2] = p;
            }
        }
        if (error)
            return false;
        for (size_t index1 = 0; index1 < neighbors.size(); ++index1) {
            for (size_t index2 = 0; index2 < MESH_SIZE; ++index2) {
                if (U[parallelarray[index1][index2].point] > parallelarray[index1][index2].value) {
                    U[parallelarray[index1][index2].point] = parallelarray[index1][index2].value;
                    Point p;

                    for (int i = 0; i < DIM; i++) {
                        p[i] = data.points[parallelarray[index1][index2].point][i];
                    }

                    if (!L_in[parallelarray[index1][index2].point]) {
                        L_in[parallelarray[index1][index2].point] = true;
                        minHeap.push(parallelarray[index1][index2].point);
                    }
                }
            }
        }
    }
    std::cout << std::endl;
    return true;
}

#endif
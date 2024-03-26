//
// Created by Enrico on 19/12/2023.
//

#if FMM

#include "EikonalSolver.hpp"
#include <queue>
#include <iostream>
#include "EikonalTraits.hpp"
#include <execution>
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

            //take time value of point p
            double p = U[i];
            //find neighbors of L[i] and get the base (the DIMENSION points with the smallest value of U
            std::vector<int> neighbors;
            std::size_t start, end;
            start = data.index.at(i).start;
            end = data.index.at(i).end;
            for (std::size_t j = start; j < end; j++) {
                neighbors.push_back(data.adjacentList[j]);
            }
            for (const auto &m_element: neighbors) {
                for (auto point: data.elements_legacy[m_element]) {
                    if (point == i || L_set[point]) continue;
                    //if point in L continue
                    //solve local problem with this point as unknown and the others as base
                    std::array<Point, MESH_SIZE> base;
                    std::size_t k = 0;
                    Eigen::Matrix<double, MESH_SIZE, 1> values;
                    for (std::size_t j = 0; j < MESH_SIZE; j++) {
                        if ((data.elements_legacy[m_element])[j] == point) {
                            base[MESH_SIZE - 1] = data.points[point];
                        } else {
                            base[k] = data.points[(data.elements_legacy[m_element])[j]];
                            values[k] = U[(data.elements_legacy[m_element])[j]];
                            k++;
                        }
                    }

                    Traits::VelocityM M;
                    M << 1.0, 0.0, 0.0,
                            0.0, 1.0, 0.0,
                            0.0, 0.0, 1.0;

                    LocalSolver<MESH_SIZE> solver(M, base, values, 5000, 10e-6);

                    double sol = solver();

                    if (sol > U[point]) continue;
                    U[point] = sol;

                    Point p;

                    for (int i = 0; i < 3; i++) {
                        p[i] = data.points[point][i];
                    }
                    if (!L_in[point]) {
                        L_in[point] = true;
                        minHeap.push(point);
                    }

                }
            }


        }
        std::cout << std::endl;
        return true;
    }
}
#endif
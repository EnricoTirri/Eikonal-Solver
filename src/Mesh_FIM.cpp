//
// Created by giorgio on 26/03/2024.
//


#if FIM

#include "EikonalSolver.hpp"
#include <queue>
#include <iostream>
#include "EikonalTraits.hpp"
#include <execution>
#include "EikonalHeapComparator.cpp"
#include "LocalSolver.hpp"

namespace Eikonal {
    using Point = typename Eikonal::Traits::Point;
#define tol 10e-6

    template<int MESH_SIZE>
    double Update(int v, const Mesh<MESH_SIZE> &data, const std::vector<double> &U) {

        double min = std::min(U[v], EikonalSolver<MESH_SIZE>::MAXF);
        //find neighbors of L[i]
        std::vector<int> neighbors;
        std::size_t start, end;
        start = data.index.at(v);
        end = data.index.at(v+1);
        for (std::size_t j = start; j < end; j++) {
            neighbors.push_back(data.adjacentList[j]);
        }


        for (auto &m_element: neighbors) {
            std::array<Point, MESH_SIZE> base;
            std::size_t k = 0;
            Eigen::Matrix<double, MESH_SIZE, 1> values;
            for (std::size_t j = 0; j < MESH_SIZE; j++) {
                if ((data.elements_legacy[m_element])[j] == v) {
                    base[MESH_SIZE - 1] = data.points[v];
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
            if (sol > U[v]) continue;
            min = std::min(min, sol);
        }
        return min;
    }

    template<int MESH_SIZE>
    void EikonalSolver<MESH_SIZE>::print_spec() {
        std::cout << "Fast Iterative Method Eikonal solver" << std::endl;
    }

    template<int MESH_SIZE>
    bool EikonalSolver<MESH_SIZE>::solve(std::vector<double> &U, const std::vector<int> &X,
                                         const Mesh<MESH_SIZE> &data) {

        for (auto &point: X) {
            if (point >= data.index.size()) {
                printf("error on initial point: %f %f does not belong to mesh\n", data.points[point].x(),
                       data.points[point].y());
                return false;
            }
        }
        std::fill(std::execution::seq, U.begin(), U.end(), MAXF);


        U.reserve(data.index.size());

        std::vector<int> L;
        std::vector<bool> L_in;
        L_in.resize(data.index.size());
        std::fill(L_in.begin(), L_in.end(), false);
        //X is the seed vertices
        //L is the active vertices
        //for every 1-ring of data in X, add the vertices to L
        for (int i : X) {
            U[i] = 0;
            //find neighbors of X[i]
            std::vector<int> neighbors;
            std::size_t start, end;
            start = data.index.at(i);
            end = data.index.at(i+1);
            for (std::size_t j = start; j < end; j++) {
                neighbors.push_back(data.adjacentList[j]);
            }
            for (auto &m_element: neighbors) {
                for (const auto &point: data.elements_legacy[m_element]) {
                    L.push_back(point);
                    L_in[point]=true;
                }
            }
        }


        std::vector<int> L_new;
        while (!L.empty()) {
            L_new.clear();
#pragma omp parallel for default(none) shared(L,L_in,L_new,U,data) num_threads(N_THREADS)
            for (auto v: L) {
                if (std::abs(U[v] - Update<MESH_SIZE>(v, data, U)) < tol) {
                    L_in[v]=false;
                    //find neighbors of L[i]
                    std::vector<int> neighbors;
                    std::size_t start, end;
                    start = data.index.at(v);
                    end = data.index.at(v+1);
                    for (std::size_t j = start; j < end; j++) {
                        neighbors.push_back(data.adjacentList[j]);
                    }

                    for (auto &m_element: neighbors) {
                        for (const auto &point: data.elements_legacy[m_element]) {
                            if (point == v || L_in[point]) continue;
                            double q = Update<MESH_SIZE>(point, data, U);
                            {
                                double p = U[point];
                                if (p > q) {
                                    U[point] = q;
#pragma omp critical
                                    L_new.push_back(point);
                                    L_in[point] = true;
                                }
                            }
                        }
                    }
                } else {
                    U[v] = Update<MESH_SIZE>(v, data, U);
#pragma omp critical
                    L_new.push_back(v);
                }
            }




        std::swap(L, L_new);
        }


        return true;
    }


}
#endif

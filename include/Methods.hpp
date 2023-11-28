//
// Created by giorgio on 28/11/23.
//

#ifndef EIKONEL_TEST_FIM_HPP
#define EIKONEL_TEST_FIM_HPP

#include "Mesh.h"
#include "SimplexData.hpp"
#include "solveEikonalLocalProblem.hpp"
#include <unordered_map>
#include <vector>

namespace methods {
    template<std::size_t DIM, std::size_t MESH_SIZE>
    struct FIM {
        typedef typename Eikonal::Eikonal_traits<DIM>::Point Point;

        bool operator()(std::unordered_map<Point, double> &U, std::vector<Point> X, std::vector<Point> L,
                        const Mesh<DIM, MESH_SIZE> &data) {

            for (auto &point: X) {
                if (!data.index.contains(point)) {
                    printf("error on initial point: %f %f does not belong to mesh\n", point.x(), point.y());
                    return false;
                }
            }

            U.reserve(data.index.size());

//#pragma omp parallel for
            for (const auto &i: data.index) {
                U.insert({i.first, MAXFLOAT});
            }
            for (const auto &i: X) {
                U[i] = 0;
                L.push_back(i);
            }

            //2. Update points in L
            while (!L.empty()) {
                std::vector<Point> next_L;
                next_L.clear();
//#pragma omp parallel for//matbe task would be better
                for (const Point &i: L) {
                    //take time value of point p
                    double p = U[i];
                    //find neighbors of L[i] and get the base (the DIMENSION points with the smallest value of U
                    std::vector<MeshElement<DIM, MESH_SIZE> *> neighbors;
                    std::size_t start, end;
                    start = data.index.at(i).start;
                    end = data.index.at(i).end;
                    for (std::size_t j = start; j < end; j++) {
                        neighbors.push_back(data.adjacentList[j]);
                    }

                    //with task maybe we can parallelize this
                    for (const auto &m_element: neighbors) {
                        for (const Point &point: *m_element) {
                            if (point == i) continue;
                            //if point in L continue
                            //solve local problem with this point as unknown and the others as base
                            std::array<Point, DIMENSION + 1> base;

                            std::size_t k = 0;
                            Eikonal::Eikonal_traits<DIMENSION>::VectorExt values;
                            for (std::size_t j = 0; j < DIMENSION + 1; j++) {
                                if ((*m_element)[j] == point) {
                                    base[DIMENSION] = point;
                                } else {
                                    base[k] = (*m_element)[j];
                                    k++;
                                }
                            }
                            for (int iter = 0; iter < DIMENSION; iter++) {
                                values[iter] = U[base[iter]];
                            }

                            Eikonal::Eikonal_traits<DIMENSION>::MMatrix M;
                            M << 1.0, 0.0,
                                    0.0, 1.0;
                            Eikonal::SimplexData<DIM> simplex{base, M};
                            Eikonal::solveEikonalLocalProblem<DIM> solver{std::move(simplex),
                                                                          values};
                            auto sol = solver();
                            //if no descent direction or no convergence kill the process
                            if (sol.status != 0) {
                                printf("error on convergence\n");
                                return false;
                            }
                            auto newU = sol.value;
/*
                    double b=-1;
                    double c=-1;
                    for(Point i2 : base)
                    {
                        if(i2!=point)
                        {
                            if (c<0)
                                c=U[i2];
                            else
                                b=U[i2];
                        }
                    }
                    auto newU2 = solveQuadratic(99,b,c);
                    if (abs(newU2-newU)>0.1)
                    {
                       // printf("difference local %f, quad: %f old: %f\n",newU,newU2,U[point]);
                    }
                    newU = newU2;*/
                            if (newU < U[point]) {
                                U[point] = newU;
                                // #pragma omp atomic
                                next_L.emplace_back(point.x(), point.y());
                            }
                        }
                    }
                }
                //remove duplicates from next_L
                std::sort(next_L.begin(), next_L.end(), [](const Point &a, const Point &b) {
                    return a[0] * 1000 + a[1] < b[0] * 1000 + b[1];
                });
                next_L.erase(std::unique(next_L.begin(), next_L.end()), next_L.end());

                L = std::move(next_L);
            }
            return true;
        }

    };

    template<std::size_t DIM, std::size_t MESH_SIZE>
    //!Fast sweeping method
    struct FSM {
        typedef typename Eikonal::Eikonal_traits<DIM>::Point Point;

        bool operator()(std::unordered_map<Point, double> &U, std::vector<Point> X, std::vector<Point> L,
                        const Mesh<DIM, MESH_SIZE> &data) {


            for (auto &point: X) {
                if (!data.index.contains(point)) {
                    printf("error on initial point: %f %f does not belong to mesh\n", point.x(), point.y());
                    return false;
                }
            }

            U.reserve(data.index.size());

//#pragma omp parallel for
            for (const auto &i: data.index) {
                U.insert({i.first, MAXFLOAT});
            }
            for (const auto &i: X) {
                U[i] = 0;
                L.push_back(i);
            }
            //heapify L, so that the smallest value is at the top
            std::make_heap(L.begin(), L.end(), [&U](const Point &a, const Point &b) {
                return U[a] < U[b];
            });


            //2. Update points in L
            while (!L.empty()) {
                Point i = L.back();
                L.pop_back();

                //take time value of point p
                double p = U[i];
                //find neighbors of L[i] and get the base (the DIMENSION points with the smallest value of U
                std::vector<MeshElement<DIM, MESH_SIZE> *> neighbors;
                std::size_t start, end;
                start = data.index.at(i).start;
                end = data.index.at(i).end;
                for (std::size_t j = start; j < end; j++) {
                    neighbors.push_back(data.adjacentList[j]);
                }

                //with task maybe we can parallelize this
                for (const auto &m_element: neighbors) {
                    for (const Point &point: *m_element) {
                        if (point == i || U[point] != MAXFLOAT) continue;
                        //if point in L continue
                        //solve local problem with this point as unknown and the others as base
                        std::array<Point, DIMENSION + 1> base;

                        std::size_t k = 0;
                        Eikonal::Eikonal_traits<DIMENSION>::VectorExt values;
                        for (std::size_t j = 0; j < DIMENSION + 1; j++) {
                            if ((*m_element)[j] == point) {
                                base[DIMENSION] = point;
                            } else {
                                base[k] = (*m_element)[j];
                                k++;
                            }
                        }
                        for (int iter = 0; iter < DIMENSION; iter++) {
                            values[iter] = U[base[iter]];
                        }

                        Eikonal::Eikonal_traits<DIMENSION>::MMatrix M;
                        M << 1.0, 0.0,
                                0.0, 1.0;
                        Eikonal::SimplexData<DIM> simplex{base, M};
                        Eikonal::solveEikonalLocalProblem<DIM> solver{std::move(simplex),
                                                                      values};
                        auto sol = solver();
                        //if no descent direction or no convergence kill the process
                        if (sol.status != 0) {
                            printf("error on convergence\n");
                            return false;
                        }
                        auto newU = sol.value;

                        if (newU < U[point]) {
                            U[point] = newU;
                            // #pragma omp atomic
                            L.emplace_back(point.x(), point.y());
                        }
                    }
                }

                //heapify L, so that the smallest value is at the top
                std::make_heap(L.begin(), L.end(), [&U](const Point &a, const Point &b) {
                    return U[a] < U[b];
                });
            }
            return true;

        }


    };


}
#endif //EIKONEL_TEST_FIM_HPP

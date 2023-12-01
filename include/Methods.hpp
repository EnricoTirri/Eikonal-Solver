//
// Created by giorgio on 28/11/23.
//

#ifndef EIKONEL_TEST_FIM_HPP
#define EIKONEL_TEST_FIM_HPP
#ifndef _NUM_THREADS
#define _NUM_THREADS 8
#endif

#include "Mesh.h"
#include "SimplexData.hpp"
#include "solveEikonalLocalProblem.hpp"
#include <Eikonal_traits.hpp>
#include <unordered_map>
#include <vector>
//#define PHDIM 3
//#define MESH_SIZE 4
#define MAXF 90000
namespace methods {
    template<int DIM, int MESH_SIZE>
    static bool FIM(std::unordered_map<typename Eikonal_traits<DIM>::Point, double> &U,
                    std::vector<typename Eikonal_traits<DIM>::Point> X,
                    std::vector<typename Eikonal_traits<DIM>::Point> L,
                    const Mesh<DIM, MESH_SIZE> &data) {
        typedef typename Eikonal_traits<DIM>::Point Point;
        for (auto &point: X) {
            if (!data.index.contains(point)) {
                printf("error on initial point: %f %f does not belong to mesh\n", point.x(), point.y());
                return false;
            }
        }

        U.reserve(data.index.size());

        for (const auto &i: data.index) {
            U[i.first] = MAXF;
        }
        for (const auto &i: X) {
            U[i] = 0;
            L.push_back(i);
        }

        //2. Update points in L
        while (!L.empty()) {
            std::vector<Point> next_L;
            next_L.clear();
            bool error = false;
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
                        std::array<Point, MESH_SIZE> base;
                        std::size_t k = 0;
                        Eigen::Matrix<double, MESH_SIZE - 1, 1> values;
                        for (int j = 0; j < DIM; ++j) {
                            values[j] = 0;
                        }
                        for (std::size_t j = 0; j < MESH_SIZE; j++) {
                            if ((*m_element)[j] == point) {
                                base[MESH_SIZE - 1] = point;
                            } else {
                                base[k] = (*m_element)[j];
                                k++;
                            }
                        }
                        for (int iter = 0; iter < DIM; iter++) {
                            values[iter] = U[base[iter]];
                        }

                        typename Eikonal::Eikonal_traits<DIM>::MMatrix M;
                        if (DIM == 2)
                            M << 1.0, 0.0,
                                    0.0, 1.0;
                        else if (DIM == 3)
                            M << 1.0, 0.0, 0.0,
                                    0.0, 1.0, 0.0,
                                    0.0, 0.0, 1.0;
                        Eikonal::SimplexData<DIM, MESH_SIZE> simplex{base, M};
                        Eikonal::solveEikonalLocalProblem<DIM, MESH_SIZE> solver{std::move(simplex),
                                                                                 values};

                        auto sol = solver();
                        //if no descent direction or no convergence kill the process
                        if (sol.status != 0) {
                            printf("error on convergence\n");
                            error = true;
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
                            Point p;
#pragma unroll
                            for (int i = 0; i < DIM; i++) {
                                p[i] = point[i];
                            }
                            next_L.emplace_back(p);
                        }
                    }
                }
            }
            if (error) return false;
            //remove duplicates from next_L
            std::sort(next_L.begin(), next_L.end(),
                      [](const Point &a, const Point &b) {
                          if constexpr (DIM == 2)
                              return a.x() * 1000000 + a.y() < b.x() * 1000000 + b.y();
                          else
                              return a.x() * 1000000 + a.y() * 1000 + a.z() < b.x() * 1000000 + b.y() * 1000 + b.z();
                      });
            next_L.erase(std::unique(next_L.begin(), next_L.end()), next_L.end());

            L = std::move(next_L);
        }
        return true;
    }


    template<int DIM, int MESH_SIZE>
    static bool FIMpp(std::unordered_map<typename Eikonal_traits<DIM>::Point, double> &U,
                      std::vector<typename Eikonal_traits<DIM>::Point> X,
                      std::vector<typename Eikonal_traits<DIM>::Point> L,
                      const Mesh<DIM, MESH_SIZE> &data) {
        typedef typename Eikonal_traits<DIM>::Point Point;
        for (auto &point: X) {
            if (!data.index.contains(point)) {
                printf("error on initial point: %f %f does not belong to mesh\n", point.x(), point.y());
                return false;
            }
        }

        U.reserve(data.index.size());

        omp_set_num_threads(_NUM_THREADS);//TODO set macro for
        for (const auto &i: data.index) {
            U[i.first] = MAXF;
        }
        for (const auto &i: X) {
            U[i] = 0;
            L.push_back(i);
        }

        //2. Update points in L
        while (!L.empty()) {
            std::vector<Point> next_L;
            next_L.clear();
            bool error = false;
#pragma omp  parallel for
            for (const Point &i: L) {
                //take time value of point p
                double p = U.at(i);
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
                        std::array<Point, MESH_SIZE> base;
                        std::size_t k = 0;
                        using VectorExt = Eigen::Matrix<double,
                                MESH_SIZE - 1, 1>;//TODO check if the dimensions are correct;

                        VectorExt values;
                        for (std::size_t j = 0; j < MESH_SIZE; j++) {
                            if ((*m_element)[j] == point) {
                                base[MESH_SIZE - 1] = point;
                            } else {
                                base[k] = (*m_element)[j];
                                k++;
                            }
                        }
                        for (int iter = 0; iter < DIM; iter++) {
                            values[iter] = U.at(base[iter]);
                        }

                        typename Eikonal::Eikonal_traits<DIM>::MMatrix M;
                        if constexpr (DIM == 2)
                            M << 1.0, 0.0,
                                    0.0, 1.0;
                        else if constexpr (DIM == 3)
                            M << 1.0, 0.0, 0.0,
                                    0.0, 1.0, 0.0,
                                    0.0, 0.0, 1.0;
                        Eikonal::SimplexData<DIM, MESH_SIZE> simplex{base, M};

                        Eikonal::solveEikonalLocalProblem<DIM, MESH_SIZE> solver{std::move(simplex),
                                                                                 values};
                        auto sol = solver();
                        //if no descent direction or no convergence kill the process
                        if (sol.status != 0) {
                            printf("error on convergence\n");
                            error = true;
                        }

#pragma omp critical (updateU)
                        {

                            auto newU = sol.value;
                            if (newU < U.at(point)) {

                                U.at(point) = newU;
                                // #pragma omp atomic
                                Point p;
#pragma omp simd
#pragma unroll
                                for (int i = 0; i < DIM; i++) {
                                    p[i] = point[i];
                                }
//#pragma omp critical (updateL)
                                next_L.emplace_back(p);
                            }
                        }
                    }
                }
            }
            if (error) return false;
            //remove duplicates from next_L
            std::sort(next_L.begin(), next_L.end(),
                      [](const Point &a, const Point &b) {
                          if (DIM == 2)
                              return a.x() * 1000000 + a.y() < b.x() * 1000000 + b.y();
                          else
                              return a.x() * 1000000 + a.y() * 1000 + a.z() < b.x() * 1000000 + b.y() * 1000 + b.z();
                      });
            next_L.erase(std::unique(next_L.begin(), next_L.end()), next_L.end());

            L = std::move(next_L);
        }
        return true;
    }



    //!Fast sweeping method
#undef DIM
#undef MESH_SIZE

    template<std::size_t DIM, std::size_t MESH_SIZE>
    bool FSM(std::unordered_map<typename Eikonal::Eikonal_traits<DIM>::Point, double> &U,
             std::vector<typename Eikonal::Eikonal_traits<DIM>::Point> X,
             std::vector<typename Eikonal::Eikonal_traits<DIM>::Point> L,
             const Mesh<DIM, MESH_SIZE> &data) {
        typedef typename Eikonal::Eikonal_traits<DIM>::Point Point;

        for (auto &point: X) {
            if (!data.index.contains(point)) {
                printf("error on initial point: %f %f does not belong to mesh\n", point.x(), point.y());
                return false;
            }
        }

        U.reserve(data.index.size());

//#pragma omp parallel for
        for (const auto &i: data.index) {
            U.insert({i.first, MAXF});
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
                    if (point == i) continue;
                    //if point in L continue
                    //solve local problem with this point as unknown and the others as base
                    std::array<Point, MESH_SIZE> base;

                    std::size_t k = 0;
                    typename Eikonal::Eikonal_traits<DIM>::VectorExt values;
                    for (std::size_t j = 0; j < MESH_SIZE; j++) {
                        if ((*m_element)[j] == point) {
                            base[MESH_SIZE - 1] = point;
                        } else {
                            base[k] = (*m_element)[j];
                            k++;
                        }
                    }
                    for (int iter = 0; iter < DIM; iter++) {
                        values[iter] = U[base[iter]];
                    }

                    typename Eikonal::Eikonal_traits<DIM>::MMatrix M;
                    if constexpr (DIM == 2)
                        M << 1.0, 0.0,
                                0.0, 1.0;
                    else if constexpr (DIM == 3)
                        M << 1.0, 0.0, 0.0,
                                0.0, 1.0, 0.0,
                                0.0, 0.0, 1.0;

                    Eikonal::SimplexData<DIM, MESH_SIZE> simplex{base, M};
                    Eikonal::solveEikonalLocalProblem<DIM, MESH_SIZE> solver{std::move(simplex),
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
                        Point p;
#pragma unroll
                        for (int i = 0; i < DIM; i++) {
                            p[i] = point[i];
                        }
                        //add element to heap L
                        L.emplace_back(p);
                        //L.emplace_back(p);
                    }
                }
            }

            //heapify L, so that the smallest value is at the top
            std::make_heap(L.begin(), L.end(), [&U](const Point &a, const Point &b) {
                return U[a] > U[b];
            });
        }
        return true;

    }


}
#endif //EIKONEL_TEST_FIM_HPP

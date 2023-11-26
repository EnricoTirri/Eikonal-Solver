#include "LocalProblem/include/SimplexData.hpp"
#include "LocalProblem/include/solveEikonalLocalProblem.hpp"
#include <iostream>
#include <algorithm>
#include <VtkParser.hpp>
#include <vector>
#include <Eigen/Core>
#include <Mesh.h>
#include <MeshLoader.hpp>

#ifndef TYPE
#define TYPE 2
#endif

#define MAXFLOAT 900000
//now we will implement this algorithm Fast iterative method (X,L)
//define hash function for Point


/*double solveQuadratic(double a, double b, double c){
    double u = c + 1;
    if (u <= b) return u;
    u = (b+c+std::sqrt(-(b*b) - (c*c) + 2*b*c +2))/2;
    if(u <= a) return u;
    u = (2*(a+b+c) + std::sqrt(4*(a+b+c)*(a+b+c) - 12*(a*a + b*b + c*c -1)))/6;
    return u;
}*/

//X will be the starting point from witch the wave will propagate, L is the active list
void FIM(std::unordered_map<Point, double> &U, std::vector<Point> X, std::vector<Point> L, const Mesh &data, FILE *fp) {

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
            std::vector<MeshElement *> neighbors;
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
                    Eikonal::SimplexData<DIMENSION> simplex{base, M};
                    Eikonal::solveEikonalLocalProblem<DIMENSION> solver{std::move(simplex),
                                                                        values};
                    auto sol = solver();
                    //if no descent direction or no convergence kill the process
                    if (sol.status != 0) {
                        printf("error on convergence\n");
                        return;
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

        fprintf(fp, "end queue\n");
        L = std::move(next_L);
    }
}


int main() {
    Mesh mesh;
    std::unordered_map<Point, double> pointData;
    std::unordered_map<MeshElement, double> elementData;

    VtkParser parser;
    parser.open("testmesh.vtk");

    if (MeshLoader::load(mesh, parser, pointData, elementData) != 1) {
        printf("unable to load the mesh");
        return 1;
    }

    printf("end parsing\n");

    //let's try it
    std::unordered_map<Point, double> U;
    std::vector<Point> L;
    std::vector<Point> X;

    X.emplace_back(0, 10);
    X.emplace_back(0, 0);

    FILE *fp = fopen("outL.txt", "w");

    time_t start = clock();
    FIM(U, X, L, mesh, fp);
    time_t end = clock();

    printf("end FIM, time elapsed: %f\n", (double) (end - start) / CLOCKS_PER_SEC);


    MeshLoader::dump(mesh, parser, U, elementData);
    parser.save("final_out.vtk");
    return 0;
}

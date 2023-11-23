#include "LocalProblem/include/SimplexData.hpp"
#include "LocalProblem/include/solveEikonalLocalProblem.hpp"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>

//subset simbol: ⊂

/*
  An eikonal equation is a non-linear first-order partial differential equation that is encountered in problems of wave propagation or in Hamiltonian system. It may be used to computer the continuos shortest path (geodesic) between points, electromagnetic potential, the arrival time of an acoustic wave, etc.
In most generic term the problem is: find the function u that satisfies an equation of the type

H(x,∇u(x))=1 , x∈Ω
u(x)=g(x) , x∈Γ⊂∂Ω

In most cases we have
H(x,∇u(x))=|∇u(x)|_M=sqrt(∇u(x)^T M(x) ∇u(x)),

where M is a symmetric positive definite function. In the simplest cases , so the problem reads simply
 |∇u|=1/c, x∈Ω
 u(x)=g(x) , x∈Γ⊂∂Ω
and c represents the celerity of the wave (often taken equal to 1). If we exclude techniques based on a regularization of the problem, we may consider two main classes of algorithms operating on a mesh covering the domain
    1. Fast marching methods. Developed originally by J. Sethian, is a special case of a level-set method. An important feature is that at each step of the algorithm the solution is updated considering among a set of “active nodes” the one which has the smallest value of u. In this way the causality principle is always satisfied (if we interpret u as the arrival time of a wave, the causality principle corresponds to the fact that the past cannot be influenced by the future).
    2. Fast sweeping methods. They can be thought as a sort of non linear Gauss-Siedel iteration, where the solution is evolved in an iterative fashion until the changes between two successive iterate is smaller than a given tolerance, or other criterions are met.
*/
/*^^


In most cases we have



where M is a symmetric positive definite function. In the simplest cases , so the problem reads simply


and c represents the celerity of the wave (often taken equal to 1). If we exclude techniques based on a regularization of the problem, we may consider two main classes of algorithms operating on a mesh covering the domain
    1. Fast marching methods. Developed originally by J. Sethian, is a special case of a level-set method. An important feature is that at each step of the algorithm the solution is updated considering among a set of “active nodes” the one which has the smallest value of u. In this way the causality principle is always satisfied (if we interpret u as the arrival time of a wave, the causality principle corresponds to the fact that the past cannot be influenced by the future).
    2. Fast sweeping methods. They can be thought as a sort of non linear Gauss-Siedel iteration, where the solution is evolved in an iterative fashion until the changes between two successive iterate is smaller than a given tolerance, or other criterions are met.
    3.
Fast marching methods have the advantage of having the solution computed after a number of steps equal to the mesh nodes, but are more difficult to parallelise than Fast sweeping type methods since the each update requires to find the active node with smallest value. They have both been developed originally for regular Cartesian meshes, but later extended to triangular and tetrahedral grids.

Both class of methods rely on the solution of a local problem, which is an optimization problem on a single mesh element. In the given literature the local problem is solved exactly, but the solution requires a sufficiently regular grid, and needs to solve a non-linear problem which is (particularly in 3D) rather nasty. For this hands-on, I suggest to solve the local problem with an optimization algorithm using the code contained in LocalProblem/,  based on modified version of the LineSearch Example of the course. The file solveEikonalProblem.cpp contains the function to use, and main_eikonal.cpp contains an example.

 */
/* What to do?
 * - Implement the algorithm described in the two references by Z. Fu et al. indicated below. Start with a scalar version. I will suggest to start with triangular meshes on a plane and tetrahedra (the algorithm si basically the same). Let alone the case of triangulated surfaces for the start. This way you can use for the Local Problem directly the solver contained in the LocalProblem/ folder.
- Move to the parallel version using OpenMP
- Use some triangular/tetrahedral mesh to test the problem. You may use the examples indicated in the generateMesh function of matlab or other tools available on the web (for instance the code triangle for triangular mesh and the code gmsh or tetgen for tetrahedral meshes). For the test you can just choose a portion of the domain boundary where to fix u and then see what happens. You may use paraview to see the solution (particularly useful in 3D)
  */

#include <vector>
#include <Eigen/Core>
#include <PointsEdge.h>
#include <memory>
#include "SimplexData.hpp"

typedef Eigen::Matrix<double, DIMENSION, 1> Point;
//now we will implement this algorithm Fast iterative method (X,L)
//define hash function for Point


double solveQuadratic(double a, double b, double c){
    double u = c + 1;
    if (u <= b) return u;
    u = (b+c+std::sqrt(-(b*b) - (c*c) + 2*b*c +2))/2;
    if(u <= a) return u;
    u = (2*(a+b+c) + std::sqrt(4*(a+b+c)*(a+b+c) - 12*(a*a + b*b + c*c -1)))/6;
    return u;
}

//X will be the starting point from witch the wave will propagate, L is the active list
void FIM(std::unordered_map<Point, double> &U, std::vector<Point> X, std::vector<Point> L, PointsEdge &data) {

    U.reserve(data.index.size());

#pragma omp parallel for
    for (auto &i: data.index) {
        U.insert({i.first, 99999e10});
    }
#pragma omp parallel for
    for (const auto &i: X) {
        U[i] = 0;
        std::size_t start = data.index[i].start;
        std::size_t end = data.index[i].end;
#pragma omp SIMD
        for (std::size_t j = start; j < end; j++) {
            L.push_back(data.adjacentList[j]);
        }

    }


    //2. Update points in L
    while (!L.empty()) {
        std::vector<Point> next_L;
#pragma omp parallel for
        for (const auto &i: L) {
            //take time value of point p
            double p = U[i];
            //find neighbors of L[i] and get the base (the DIMENSION points with the smallest value of U)
            std::vector<Point> neighbors;
            std::size_t start, end;
            start = data.index[i].start;
            end = data.index[i].end;
            for (std::size_t j = start; j < end; j++) {
                neighbors.push_back(data.adjacentList[j]);
            }
            //sorting neighbors on arrival time
            std::sort(neighbors.begin(), neighbors.end(), [&U](Point const &a, Point const &b) {
                return U[a] < U[b];
            });
            //build base
            std::array<Point, DIMENSION + 1> base;
            std::size_t j;
            Point last;
            for (j = 0; j < DIMENSION && j < (neighbors.size() - 1); j++) {
                base[j] = neighbors[j];
                last = neighbors[j];
            }
            //if j<DIMENSION-1 the last points will be duplicates
            for (std::size_t k = j; k < DIMENSION; k++) {
                base[k] = last;
            }
            base[DIMENSION] = i;

            //test M = identity
            Eikonal::Eikonal_traits<DIMENSION>::MMatrix M;
            //TODO we need to understand how to construct a useful M matrix with the speeds
            M = Eikonal::Eikonal_traits<DIMENSION>::MMatrix::Identity();

            //build data structure for localproblemsolver
            Eikonal::SimplexData<DIMENSION> simplex{base, M};

            using VectorExt = Eikonal::Eikonal_traits<DIMENSION>::VectorExt;
            VectorExt values;
            //init values size to DIMENSION
            values.resize(DIMENSION);
            //value will be the value of U at the base points - the last one
            int j2;
            double last_v;
            for (j2 = 0; j2 < j; j2++) {
                values[j2] = U[base[j2]];
                last_v = values[j2];
            }
            for (std::size_t k = j2; k < DIMENSION; k++) {
                values[k] = last_v;
            }


            Eikonal::solveEikonalLocalProblem<DIMENSION> solver{std::move(simplex),
                                                            values};//TODO: edit local solver to ignore value negative or to allow infty
            auto sol = solver();
            //if no descent direction or no convergence kill the process
            if (sol.status != 0) {
                printf("error on convergence");
                return;
            }
            auto newU = sol.value;

            //newU = solveQuadratic(999999,values[1],values[0]);
            //  printf("%f\n",newU);
            //if the new value is smaller than the old one update it
            if (newU < p) {
                U[i] = newU;
                /*TODO: for each neighbor of L[i] check if the propagation would improve time, if so add it to L
                 * maybe we don't need to check because we are already doing it in the while loop ? not sure
                 * FOR NOW IT SEAMS TO WORK
                 */

                for (const auto &k: neighbors) {
#pragma omp atomic
                    next_L.push_back(k);
                }
            }

        }

        L = std::move(next_L);
    }
}


int main() {
    //open file
    FILE *fp;
    fp = fopen("test1.txt", "r");
    if (fp == NULL) {
        printf("Test file not found");
        return 1;
    }

    //read # points
    int n_points;
    fscanf(fp, "%d", &n_points);
    //read # edges
    int n_edges;
    fscanf(fp, "%d", &n_edges);

    //support structure for parsing file
    std::unordered_map<Point, std::vector<Point>> readed;

    //final structure
    PointsEdge pointsEdge;

    //parsing function
    while (!feof(fp)) {
        Point p1, p2;

        (void) fscanf(fp, "%lf,%lf %lf,%lf", &p1[0], &p1[1], &p2[0], &p2[1]);
        readed[p1].push_back(p2);
        readed[p2].push_back(p1);
    }
    printf("end parsing");

    //remove duplicates in each vector in map
//    for (auto &i: readed) {
//        std::sort(i.second.begin(), i.second.end(), [](Point const &a, Point const &b) {
//            return (pow(a[0], 2) + pow(a[1], 2)) < (pow(b[0], 2) + pow(b[1], 2));
//        });
//        i.second.erase(std::unique(i.second.begin(), i.second.end()), i.second.end());
//    }
    int k = 0;
    size_t prev = 0;
    for (auto &i: readed) {
        if (k == 0)
            pointsEdge.index[i.first].start = 0;
        else
            pointsEdge.index[i.first].start = prev;

        for (const auto &j: i.second) {
            pointsEdge.adjacentList.push_back(j);
            k++;
        }
        pointsEdge.index[i.first].end = k;
        prev = k;
    }
    printf("t");

    printf("t");

    //let's try it
    std::unordered_map<Point, double> U;
    std::vector<Point> L;
    std::vector<Point> X;
    Point start;
    start << 0.0, 0.0;
    X.push_back(start);
    FIM(U, X, L, pointsEdge);

    double square[16][16];
    std::ofstream out("out.csv");

    out << "x" << "," << "y" << "," << "z" << std::endl;
    for (auto a: U) {
        out << a.first[0] << "," << a.first[1] << "," << a.second << std::endl;
    }
//
//    for (int i = 0; i < 16; i++) {
//        int j;
//        for (j = 0; j < 15; j++)
//            out << square[i][j] << "\t";
//        out << square[i][j] << std::endl;
//    }

    out.flush();
    out.close();

//    for(auto i: U)
//    {
//        square[(floor(i.first[0])),
//    }


    return 0;
}

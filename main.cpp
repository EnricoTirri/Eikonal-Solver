#include <iostream>

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
#define PHDIM 2

#include <vector>
#include <Eigen/Core>
#include "Eigen/Core"
#include "Eigen/Core"
#include <PointsEdge.h>
#include <memory>

typedef Eigen::Matrix<double, PHDIM, 1> Point;
//now we will implement this algorithm Fast iterative method (X,L)
//define hash function for Point


//X will be the starting point from witch the wave will propagate, L is the active list
void FIM(std::unordered_map<Point, double> &U, std::vector<Point> X, std::vector<Point> L, PointsEdge &data) {
    //1. Initialization (X : a set of grid nodes, L : active list)
    U.reserve(data.index.size());
    std::vector<Point> startPoints;

    for (auto &i: data.index) {
        //if i.first is in X then U(i.first)=0 else U(i.first)=inf
        if (std::find(X.begin(), X.end(), i.first) != X.end()) {
            U.insert({i.first, 0});
            startPoints.push_back(i.first);
        } else {
            U.insert({i.first, std::numeric_limits<double>::infinity()});
        }
    }
    std::size_t start, end;
    for (const auto &i: startPoints) {
        start = data.index[i].start;
        end = data.index[i].end;
        for (std::size_t j = start; j < end; j++) {
            L.push_back(data.adjacentList[j]);
        }

    }


    //2. Update points in L
    while (!L.empty()) {
        //for every point in L
        for (const auto &i: L) {
            double p = U[i];
            //find neighbors of L[i] and get the base (the PHDIM points with the smallest value of U)
            std::vector<Point> neighbors;
            std::size_t start, end;
            start = data.index[i].start;
            end = data.index[i].end;
            for (std::size_t j = start; j < end; j++) {
                neighbors.push_back(data.adjacentList[j]);
            }
            std::sort(neighbors.begin(), neighbors.end(), [&U](Point const &a, Point const &b) {
                return U[a] < U[b];
            });
            std::vector<Point> base;
            for (int j = 0; j < PHDIM; j++) {
                base.push_back(neighbors[j]);
            }
            //solve the local problem
            //TODO: solve the local problem with LocalProblem
            //TODO: for each neighbor of L[i] check if the propagation would improve time, if so add it to L



        }


    }
}


int main() {
    FILE *fp;
    fp = fopen("test1.txt", "r");
    int n_points;
    fscanf(fp, "%d", &n_points);
    int n_edges;
    fscanf(fp, "%d", &n_edges);
    std::unordered_map<Point, std::vector<Point>> readed;
    std::vector<Point> Points;
    std::vector<Point> Edges;
    PointsEdge pointsEdge;
    //smart pointer
    std::unique_ptr<PointsEdge> p;
    p = std::make_unique<PointsEdge>();

    while (!feof(fp)) {
        Point p1, p2;

        (void) fscanf(fp, "%lf ,%lf %lf ,%lf", &p1[0], &p1[1], &p2[0], &p2[1]);
        readed[p1].push_back(p2);
        readed[p2].push_back(p1);
    }

    printf("t");
    //remove duplicates in each vector in map
    for (auto &i: readed) {
        std::sort(i.second.begin(), i.second.end(), [](Point const &a, Point const &b) {
            return (pow(a[0], 2) + pow(a[1], 2)) < (pow(b[0], 2) + pow(b[1], 2));
        });
        i.second.erase(std::unique(i.second.begin(), i.second.end()), i.second.end());
    }
    int k = 0;
    for (auto &i: readed) {
        pointsEdge.index[i.first].start = 0;
        for (const auto &j: i.second) {
            pointsEdge.adjacentList.push_back(j);
            k++;
        }
        pointsEdge.index[i.first].end = k;
    }
    printf("t");
    for (auto &i: readed) {
        pointsEdge.list.push_back(i.first);
    }
    printf("t");
   return 0;
}

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

typedef Eigen::Matrix<double, PHDIM, 1> Point;
//now we will implement this algorithm Fast iterative method (X,L)
//define hash function for Point
namespace std {
    template<>
    struct hash<Point> {
        std::size_t operator()(const Point &k) const {
            using std::size_t;
            using std::hash;
            using std::string;
            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:
            return ((hash<double>()(k[0])
                     ^ (hash<double>()(k[1]) << 1)) >> 1);
        }
    };
}

void FIM(std::unordered_map<Point, double> &U, std::vector<Point> X, std::vector<Point> L,
         std::unordered_map<Point, std::vector<Point>> neighbors) {
    //1. Initialization (X : a set of grid nodes, L : active list)
    U.reserve(X.size());


    for (auto &i: X) {
        //if x[i] is the 0,0 point
        if (i[0] == 0 && i[1] == 0) {
            U.insert({i, 0});
        } else {
            U.insert({i, std::numeric_limits<double>::infinity()});
        }

    }
    for (int i = 0; i < X.size(); i++) {
        for (int j = 0; j < X[i].size(); j++) {
            if (X[i][j] == 0) {
                L.emplace_back(i, j);
            }
        }
    }
    //2. Update points in L
    while (!L.empty()) {
        //for every point in L
        for (const auto &i: L) {
            double p = U[i];
            //find neighbors of L[i] and get the base
            std::vector<Point> neighbors2 = neighbors[i];


        }


        for (int i = 0; i < L.size(); i++) {
            for (int j = 0; j < L[i].size(); j++) {
                double p, q;
                //double p =U
                //p=U(x)
                //q=solution of g(x)=0
                //U(x)=q
                if (abs(p - q) < 0.0001) {
                    for (int i = 0; i < X.size(); i++) {
                        for (int j = 0; j < X[i].size(); j++) {
                            if (X[i][j] == 0) {
                                //add x to L
                            }
                        }
                    }
                }
                for (int i = 0; i < X.size(); i++) {
                    for (int j = 0; j < X[i].size(); j++) {
                        if (X[i][j] == 0) {
                            //U(x)=q
                        }
                    }
                }
            }
        }
    }
}

int main() {
   std::cout<<"Hello world" <<std::endl;
   return 0;





}

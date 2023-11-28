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
#undef DIMENSION

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


int main() {
    Mesh mesh;
    std::unordered_map<Point, double> pointData;
    std::unordered_map<MeshElement, double> elementData;

    VtkParser parser;
    parser.open("smesh.vtk");

    if (MeshLoader::load(mesh, parser, pointData, elementData) != 1) {
        printf("unable to load the mesh");
        return 1;
    }

    printf("end parsing\n");

    //let's try it
    std::unordered_map<Point, double> U;
    std::vector<Point> L;
    std::vector<Point> X;

    X.emplace_back(0, 0);
    X.emplace_back(21,0);
    X.emplace_back(10.5,0);

    time_t start = clock();
    bool success = FIM(U, X, L, mesh);
    time_t end = clock();

    if(success) {
        printf("end FIM, time elapsed: %f\n", (double) (end - start) / CLOCKS_PER_SEC);

        MeshLoader::dump(mesh, parser, U, elementData);
        parser.save("final_out.vtk");
    }
    return 0;
}

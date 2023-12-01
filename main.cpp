#include "LocalProblem/include/SimplexData.hpp"
#include "LocalProblem/include/solveEikonalLocalProblem.hpp"
#include "Methods.hpp"
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
    const size_t pointdim = 3;
    const size_t meshdim = 3;
    const string filename = "testmesh";
    {
        Mesh<pointdim, meshdim> mesh;
        using Point = Eikonal_traits<pointdim>::Point;
        std::unordered_map<Point, double> pointData;
        std::unordered_map<MeshElement<pointdim, meshdim>, double> elementData;

        VtkParser parser;
        parser.open(filename + ".vtk");

        if (MeshLoader<pointdim, meshdim>::load(mesh, parser, pointData, elementData) != 1) {
            printf("unable to load the mesh");
            return 1;
        }

        printf("end parsing\n");

        //let's try it
        std::unordered_map<Point, double> U;
        std::vector<Point> L;
        std::vector<Point> X;
        //-0.0378297 0.12794 0.00447467

        // X.emplace_back(-0.0378297, 0.12794, 0.00447467);

        X.emplace_back(0, 0, 0);
        X.emplace_back(1, 1, 0);
        time_t start2 = clock();

        bool success = methods::FIM<pointdim, meshdim>(U, X, L, mesh);
        time_t end2 = clock();

        if (success) {
            printf("end FIM, time elapsed: %f\n", (double) (end2 - start2) / CLOCKS_PER_SEC);

            MeshLoader<pointdim, meshdim>::dump(mesh, parser, U, elementData);
            parser.save(filename + "_out.vtk");
        } else {
            printf("FIM failed\n");
        }
        U.clear();
        L.clear();
    }

    {
        Mesh<3u, 3u> mesh;
        using Point = Eikonal_traits<3u>::Point;
        std::unordered_map<Point, double> pointData;
        std::unordered_map<MeshElement<3, 3>, double> elementData;

        VtkParser parser;
        parser.open("bunny.vtk");

        if (MeshLoader<3, 3>::load(mesh, parser, pointData, elementData) != 1) {
            printf("unable to load the mesh");
            return 1;
        }

        printf("end parsing\n");

        //let's try it
        std::unordered_map<Point, double> U;
        std::vector<Point> L;
        std::vector<Point> X;
        //-0.0378297 0.12794 0.00447467
        X.emplace_back(-0.0378297, 0.12794, 0.00447467);

        // X.emplace_back(0, 0, 0);
        timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        bool success = methods::FIMpp<3, 3>(U, X, L, mesh);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (success) {
            double elapsed = (end.tv_sec - start.tv_sec);
            elapsed += (end.tv_nsec - start.tv_nsec) / 1000000000.0;
            printf("end FIM, time elapsed: %f\n", elapsed);

            MeshLoader<3, 3>::dump(mesh, parser, U, elementData);
            parser.save("bunny_outpp.vtk");
        } else {
            printf("FIM failed\n");
        }
        U.clear();
        L.clear();
    }
    //let's try the fast sweeping method

    return 0;
}

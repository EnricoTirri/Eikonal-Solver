#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"
#include "LocalProblem/include/SimplexData.hpp"
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


int main(int argc, char *argv[]) {
    if (argc < 6) {
        printf("Usage: %s <filename.vtk> <output.vtk> <pointdim> <meshdim> <id1> ... <idn>\n",
               argv[0]);
        return 1;
    }
    bool success = false;
    const char *filename = argv[1];
    const char *output_filename = argv[2];
    const size_t pointdim = std::atoi(argv[3]);
    const size_t meshdim = std::atoi(argv[4]);

    std::vector<int> startingPoints;
    for (std::size_t i = 5; i < argc; i += pointdim) {
        int point;
        point = std::atoi(argv[i]);
        startingPoints.push_back(point);
    }
    //try to open the file
    FILE *file = fopen(filename, "r");
    if (file == nullptr) {
        printf("unable to open the file %s\n", filename);
        return 1;
    }
    fclose(file);
    if (pointdim != 2 && pointdim != 3) {
        printf("pointdim must be 2 or 3\n");
        return 1;
    }
    if (meshdim != 3 && meshdim != 4) {
        printf("meshdim must be 3 or 4\n");
        return 1;
    }
    if (pointdim == 2 && meshdim != 3) {
        printf("meshdim must be 3 if pointdim is 2\n");
        return 1;
    }


    VtkParser parser;
    parser.open(filename);

    //correct inputs 2,3 3,3 3,4
    if (pointdim == 2) {

        Mesh<2u, 3u> mesh;
        using Point = Eikonal_traits<2u>::Point;
        std::vector<double> pointData;
        std::vector<double> elementData;

        if (MeshLoader<2, 3>::load(mesh, parser, pointData, elementData) != 1) {
            printf("unable to load the mesh");
            return 1;
        }
        parser = VtkParser();

        printf("end parsing\n");

        //let's try it
        std::vector<double> U(mesh.points.size());
        std::vector<int> X;

        for (const auto point: startingPoints) {


            X.emplace_back(point);

        }

        timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        success = methods::FIM<2, 3>(U, X, mesh);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (success) {
            auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
            elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
            printf("end FIM, time elapsed: %f\n", elapsed);


            MeshLoader<2, 3>::dump(mesh, parser, U, elementData);
            parser.save(output_filename);
        } else {
            printf("FIM failed\n");
        }
        U.clear();
    } else if (meshdim == 4)//3,4
    {
        Mesh<3u, 4u> mesh;
        using Point = Eikonal_traits<3u>::Point;
        std::vector<double> pointData;
        std::vector<double> elementData;

        if (MeshLoader<3, 4>::load(mesh, parser, pointData, elementData) != 1) {
            printf("unable to load the mesh");
            return 1;
        }

        parser = VtkParser();

        printf("end parsing\n");

        //let's try it
        std::vector<double> U(mesh.points.size());

        std::vector<int> X;

        for (auto point: startingPoints) {

            X.emplace_back(point);
        }

        timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        success = methods::FIM<3, 4>(U, X, mesh);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (success) {
            auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
            elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
            printf("end FIM, time elapsed: %f\n", elapsed);

            parser = VtkParser();
            MeshLoader<3, 4>::dump(mesh, parser, U, elementData);
            parser.save(output_filename);
        } else {
            printf("FIM failed\n");
        }
        U.clear();
    } else //3,3
    {
        Mesh<3u, 3u> mesh;
        using Point = Eikonal_traits<3u>::Point;
        std::vector<double> pointData;
        std::vector<double> elementData;

        if (MeshLoader<3, 3>::load(mesh, parser, pointData, elementData) != 1) {
            printf("unable to load the mesh");
            return 1;
        }

        parser = VtkParser();

        printf("end parsing\n");

        //let's try it
        std::vector<double> U(mesh.points.size());
        std::vector<int> X;

        for (auto &point: startingPoints) {

            X.emplace_back(point);
        }

        timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        success = methods::FIM<3, 3>(U, X, mesh);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (success) {
            auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
            elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
            printf("end FIM, time elapsed: %f\n", elapsed);


            MeshLoader<3, 3>::dump(mesh, parser, U, elementData);
            parser.save(output_filename);
        } else {
            printf("FIM failed\n");
        }
        U.clear();
    }

    return 0;
}

#pragma clang diagnostic pop
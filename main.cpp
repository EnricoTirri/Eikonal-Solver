#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"

#include <iostream>
#include "VtkParser.hpp"
#include <vector>
#include "Mesh.hpp"
#include "MeshLoader.hpp"
#include "EikonalSolver.hpp"

using namespace Eikonal;

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);

    if (argc < 5) {
        printf("Usage: %s <filename.vtk> <output.vtk> <meshdim> <id1> [<id2> ...]\n",
               argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    const char *output_filename = argv[2];
    const size_t meshdim = std::atoi(argv[3]);

    std::vector<int> startingPoints;
    for (std::size_t i = 4; i < argc; ++i) {
        int point;
        point = std::atoi(argv[i]);
        startingPoints.push_back(point);
    }

    // try to open the file
    FILE *file = fopen(filename, "r");
    if (file == nullptr) {
        printf("unable to open the file %s\n", filename);
        return 1;
    }
    fclose(file);

    if (meshdim != 3 && meshdim != 4) {
        printf("meshdim must be 3 or 4\n");
        return 1;
    }

    VtkParser parser = VtkParser();
    parser.open(filename);

    bool success;


    std::vector<double> pointData;
    std::vector<double> elementData;
    std::vector<double> U;
    std::vector<int> X;
    auto prepare =0.0;
    auto compute = 0.0;

    if (meshdim == 4) // 3,4
    {
        Mesh<4u> mesh;

        MeshLoader<4> loader;
        std::cout << "LOADING MESH...\t";
        if (loader.load(mesh, parser, pointData, elementData) != 1) {
            std::cout << "unable to load the mesh" << std::endl;
            return 1;
        } else {
            std::cout << "loaded" << std::endl;
        }

        for (auto point: startingPoints) {
            X.emplace_back(point);
        }

        EikonalSolver<4> solver;

        std::cout << "LAUNCHING SOLVER ..." << std::endl;
        solver.print_spec();

        success = solver.solve(U, X, mesh);
        prepare = solver.prepare;
        compute = solver.compute;

        loader.dump(mesh, parser, U, elementData);
        U.clear();

    } else // 3,3
    {
        Mesh<3u> mesh;

        MeshLoader<3> loader;
        std::cout << "LOADING MESH...\t";
        if (loader.load(mesh, parser, pointData, elementData) != 1) {
            std::cout << "unable to load the mesh" << std::endl;
            return 1;
        } else {
            std::cout << "loaded" << std::endl;
        }

        for (auto &point: startingPoints) {
            X.emplace_back(point);
        }

        EikonalSolver<3> solver;
        std::cout << "LAUNCHING SOLVER ..." << std::endl;
        solver.print_spec();


        success = solver.solve(U, X, mesh);
        prepare = solver.prepare;
        compute = solver.compute;

        loader.dump(mesh, parser, U, elementData);
        U.clear();
    }

    if (success) {
        std::cout << "solver succeeded, prepare time: " << prepare/ 1000000000.0 << " compute time: " << compute/1000000000.0 << std::endl;

        parser.save(output_filename);
    } else {
        std::cout << "solver failed" << std::endl;
    }

    return 0;
}

#pragma clang diagnostic pop
#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"

#include <iostream>
#include <VtkParser.hpp>
#include <vector>
#include <Mesh.h>
#include <MeshLoader.hpp>
#include <EikonalSolver.hpp>

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

    VtkParser parser;
    parser.open(filename);

    bool success;


    if (meshdim == 4) // 3,4
    {
        Mesh<4u> mesh;
        std::vector<double> pointData;
        std::vector<double> elementData;

        MeshLoader<4> loader;
        if (loader.load(mesh, parser, pointData, elementData) != 1) {
            printf("unable to load the mesh");
            return 1;
        }

        parser = VtkParser();


        std::vector<double> U(mesh.points.size());

        std::vector<int> X;

        for (auto point: startingPoints) {
            X.emplace_back(point);
        }

        timespec start, end;
        EikonalSolver<4> solver;

        solver.print_spec();

        clock_gettime(CLOCK_MONOTONIC, &start);
        success = solver.solve(U, X, mesh);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (success) {
            auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
            elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
            printf("end solver, time elapsed: %f\n", elapsed);

            parser = VtkParser();
            loader.dump(mesh, parser, U, elementData);
            parser.save(output_filename);
        } else {
            printf("solver failed\n");
        }
        U.clear();
    } else // 3,3
    {

        Mesh<3u> mesh;
        std::vector<double> pointData;
        std::vector<double> elementData;

        MeshLoader<3> loader;
        if (loader.load(mesh, parser, pointData, elementData) != 1) {
            printf("unable to load the mesh");
            return 1;
        }

        parser = VtkParser();


        std::vector<double> U(mesh.points.size());
        std::vector<int> X;

        for (auto &point: startingPoints) {

            X.emplace_back(point);
        }

        timespec start, end;
        EikonalSolver<3> solver;

        solver.print_spec();

        clock_gettime(CLOCK_MONOTONIC, &start);
        success = solver.solve(U, X, mesh);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (success) {
            auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
            elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
            printf("end solver, time elapsed: %f\n", elapsed);

            loader.dump(mesh, parser, U, elementData);
            parser.save(output_filename);
        } else {
            printf("solver failed\n");
        }
        U.clear();
    }

    return 0;
}

#pragma clang diagnostic pop
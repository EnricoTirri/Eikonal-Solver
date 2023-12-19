#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"

#include <iostream>
#include <VtkParser.hpp>
#include <vector>
#include <Mesh.h>
#include <MeshLoader.hpp>
#include <EikonalSolver.hpp>

#ifndef TYPE
#define TYPE 2
#endif

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);

    if (argc < 6) {
        printf("Usage: %s <filename.vtk> <output.vtk> <pointdim> <meshdim> <id1> [<id2> ...]\n",
               argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    const char *output_filename = argv[2];
    const size_t pointdim = std::atoi(argv[3]);
    const size_t meshdim = std::atoi(argv[4]);

    std::vector<int> startingPoints;
    for (std::size_t i = 5; i < argc; ++i) {
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

    bool success;

    if (pointdim == 2) {

        Mesh<2u, 3u> mesh;
        std::vector<double> pointData;
        std::vector<double> elementData;

        if (MeshLoader<2, 3>::load(mesh, parser, pointData, elementData) != 1) {
            printf("unable to load the mesh");
            return 1;
        }
        parser = VtkParser();

        std::vector<double> U(mesh.points.size());
        std::vector<int> X;

        for (const auto point: startingPoints) {
            X.emplace_back(point);
            printf("starting point: %f %f\n", mesh.points[point][0], mesh.points[point][1]);
        }

        timespec start, end;
        EikonalSolver<2, 3> solver;

        solver.print_spec();

        clock_gettime(CLOCK_MONOTONIC, &start);
        success = solver.solve(U, X, mesh);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (success) {
            auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
            elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
            printf("end solver, time elapsed: %f\n", elapsed);

            MeshLoader<2, 3>::dump(mesh, parser, U, elementData);
            parser.save(output_filename);
        } else {
            printf("solver failed\n");
        }
        U.clear();
    } else if (meshdim == 4) // 3,4
    {
        Mesh<3u, 4u> mesh;
        std::vector<double> pointData;
        std::vector<double> elementData;

        if (MeshLoader<3, 4>::load(mesh, parser, pointData, elementData) != 1) {
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
        EikonalSolver<3, 4> solver;

        solver.print_spec();

        clock_gettime(CLOCK_MONOTONIC, &start);
        success = solver.solve(U, X, mesh);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (success) {
            auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
            elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
            printf("end solver, time elapsed: %f\n", elapsed);

            parser = VtkParser();
            MeshLoader<3, 4>::dump(mesh, parser, U, elementData);
            parser.save(output_filename);
        } else {
            printf("solver failed\n");
        }
        U.clear();
    } else // 3,3
    {

        Mesh<3u, 3u> mesh;
        std::vector<double> pointData;
        std::vector<double> elementData;

        if (MeshLoader<3, 3>::load(mesh, parser, pointData, elementData) != 1) {
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
        EikonalSolver<3, 3> solver;

        solver.print_spec();

        clock_gettime(CLOCK_MONOTONIC, &start);
        success = solver.solve(U, X, mesh);
        clock_gettime(CLOCK_MONOTONIC, &end);

        if (success) {
            auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
            elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
            printf("end solver, time elapsed: %f\n", elapsed);

            MeshLoader<3, 3>::dump(mesh, parser, U, elementData);
            parser.save(output_filename);
        } else {
            printf("solver failed\n");
        }
        U.clear();
    }

    return 0;
}

#pragma clang diagnostic pop
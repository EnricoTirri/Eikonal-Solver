#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"

#include <iostream>
#include "VtkParser.hpp"
#include <vector>
#include "Mesh.hpp"
#include "MeshLoader.hpp"
#include "EikonalSolver.hpp"
#include <fstream>
#include <dirent.h>
#include <filesystem>

using namespace Eikonal;

#define ITERATIONS 10

void exportCsv(const std::string &filename, const std::vector<std::string> &results) {
    std::ofstream out(filename);
    if (out.bad()) {
        std::cout << "Bad output filename" << std::endl;
        out.open("./temp_ris_save.vtk");
    }

    out << "N points,N elements,Prepare avg time,Compute avg time,Best total time" << std::endl;
    for(auto ris : results)
        out << ris << std::endl;
}

std::string to_csv_line(const int&n_points, const int &n_elements, const double &prepare_time, const double &compute_time, const double &best_total){
   return (std::to_string(n_points) + "," + std::to_string(n_elements) + "," + std::to_string(prepare_time) + "," + std::to_string(compute_time) + ',' + std::to_string(best_total));
}

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);

    if (argc < 4) {
        printf("Usage: %s <input_directory> <output_ris.csv> <meshdim>\n",
               argv[0]);
        return 1;
    }

    const char *dir_name = argv[1];
    const char *output_filename = argv[2];
    const size_t meshdim = std::atoi(argv[3]);

    std::vector<int> startingPoints;
    startingPoints.push_back(0);
//    for (std::size_t i = 4; i < argc; ++i) {
//        int point;
//        point = std::atoi(argv[i]);
//        startingPoints.push_back(point);
//    }

    if (meshdim != 3 && meshdim != 4) {
        printf("meshdim must be 3 or 4\n");
        return 1;
    }

    std::filesystem::directory_iterator directory(dir_name);

    std::vector<std::string> results;

    for (const auto &entry: directory) {
        const char *filename = entry.path().c_str();
        // try to open the file
        FILE *file = fopen(filename, "r");
        if (file == nullptr) {
            std::cout << "UNABLE TO OPEN THE FILE : " << filename << std::endl;
            std::cout << "SKIPPING FILE ... " << std::endl;
            continue;
        }else{
            std::cout << "TESTING FILE : " << filename << std::endl;
            std::cout << "--------------------------------------" << std::endl;
        }
        fclose(file);

        VtkParser parser = VtkParser();
        parser.open(filename);

        bool global_success = true;


        std::vector<double> pointData;
        std::vector<double> elementData;
        std::vector<double> U;
        std::vector<int> X;



        double compute_media = 0;
        double prepare_media = 0;
        double best_total = MAXFLOAT;
        int n_points = 0;
        int n_elements = 0;

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
            n_points = mesh.points.size();
            n_elements = mesh.adjElementPtr.back();
            std::cout << "TOTAL POINTS : " << n_points << std::endl;
            std::cout << "TOTAL ELEMENTS : " << n_elements << std::endl;
            std::cout << "MESH TYPE : TETRAHEDRAL" << std::endl;
            std::cout << "--------------------------------------" << std::endl;
            for (int i = 0; i < ITERATIONS; ++i) {
                std::cout << "ITERATION : " << i << std::endl;
                U.clear();
                bool success = solver.solve(U, X, mesh);
                global_success &= success;
                std::cout << "PREPARE TIME : " << solver.prepare / 1000000000.0 << std::endl;
                std::cout << "COMPUTE TIME : " << solver.compute / 1000000000.0 << std::endl;
                std::string ris = success ? "SUCCEEDED" : "FAILED";
                std::cout << "SOLVER STATUS : " << ris << std::endl;
                std::cout << "--------------------------------------" << std::endl;
                compute_media += solver.compute;
                prepare_media += solver.prepare;
                if(best_total > solver.compute + solver.prepare) best_total = solver.compute + solver.prepare;
            }

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
            n_points = mesh.points.size();
            n_elements = mesh.adjElementPtr.back();
            std::cout << "TOTAL POINTS : " << n_points << std::endl;
            std::cout << "TOTAL ELEMENTS : " << n_elements << std::endl;
            std::cout << "MESH TYPE : TRIANGULAR" << std::endl;
            std::cout << "--------------------------------------" << std::endl;
            for (int i = 0; i < ITERATIONS; ++i) {
                std::cout << "ITERATION : " << i << std::endl;
                U.clear();
                bool success = solver.solve(U, X, mesh);
                global_success &= success;
                std::cout << "PREPARE TIME : " << solver.prepare / 1000000000.0 << std::endl;
                std::cout << "COMPUTE TIME : " << solver.compute / 1000000000.0 << std::endl;
                std::string ris = success ? "SUCCEEDED" : "FAILED";
                std::cout << "SOLVER STATUS : " << ris << std::endl;
                std::cout << "--------------------------------------" << std::endl;
                compute_media += solver.compute;
                prepare_media += solver.prepare;
                if(best_total > solver.compute + solver.prepare) best_total = solver.compute + solver.prepare;
            }

        }

        if (global_success) {
            std::cout << "SOLVER GLOBALLY SUCCEEDED - RESULTS :" << std::endl;
            prepare_media /= (1000000000.0 * ITERATIONS);
            compute_media /= (1000000000.0 * ITERATIONS);
            best_total /= (1000000000.0);
            std::cout << "MEDIA PREPARE TIME : " << prepare_media << std::endl;
            std::cout << "MEDIA COMPUTE TIME : " << compute_media << std::endl;
            std::cout << "BEST TOTAL TIME : " << best_total << std::endl;
            std::string res = to_csv_line(n_points, n_elements, prepare_media, compute_media, best_total);
            results.push_back(res);
        } else {
            std::cout << "SOLVER GLOBALLY FAILED" << std::endl;
        }

        std::cout << "--------------------------------------" << std::endl;
    }

    exportCsv(output_filename, results);
    return 0;
}

#pragma clang diagnostic pop
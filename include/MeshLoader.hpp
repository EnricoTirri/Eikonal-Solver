//
// Created by Enrico on 25/11/2023.
//

#ifndef EIKONEL_TEST_MESHLOADER_H
#define EIKONEL_TEST_MESHLOADER_H



#include <Mesh.h>
#include <Eikonal_traits.hpp>
#include <VtkParser.hpp>

using namespace Eikonal;
using namespace std;

template<size_t DIM, size_t MESHSIZE>
//#define DIM 3
//#define MESHSIZE 3
struct MeshLoader {

    typedef Mesh<DIM,MESHSIZE>::MeshElement M;
    typedef Eikonal_traits<DIM>::Point P;

    static int load(Mesh<DIM, MESHSIZE> &mesh, const VtkParser &parser,
                    vector<double> &pointsData, vector<double> &elementData) {
        if (parser.status != 1)
            return 0;

#ifdef MSHLOADER_VERBOSE
        timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        printf("\n(MSHLOADER)[LOAD]: Start loading from parser\n");
#endif

        mesh.index.clear();
        mesh.elements.clear();
        mesh.adjacentList.clear();

        //support structure for parsing file
        std::vector<std::vector<int>> read;
        read.resize(parser.points.size());
#ifdef MSHLOADER_VERBOSE
        cout << "[LOAD]: starting loading points ... ";
        cout.flush();
#endif


        for (auto point: parser.points) {
            P &p = mesh.points.emplace_back();
#pragma unroll DIM
            for (int i = 0; i < DIM; ++i)
                p[i] = point.vec[i];
        }

        for (auto cell: parser.cells) {
            M &m = mesh.elements.emplace_back();
#pragma unroll MESHSIZE
            for (int j = 0; j < MESHSIZE; ++j)
                m[j] = cell.point_ids[j];
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[LOAD]: starting loading mesh elements ... ";
        cout.flush();
#endif

        for (int j = 0; j < mesh.elements.size(); ++j) {
            for (int i = 0; i < MESHSIZE; i++) {
                read[mesh.elements[j][i]].push_back(j);
            }
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[LOAD]: starting loading adjacents ... ";
        cout.flush();
#endif

        size_t k = 0;
        size_t prev = 0;
        mesh.index.resize(parser.points.size());

        for (int i = 0; i < read.size(); ++i) {
            if (k == 0)
                mesh.index[i].start = 0;
            else
                mesh.index[i].start = prev;

            for (auto j: read[i]) {
                mesh.adjacentList.push_back(j);
                k++;
            }
            mesh.index[i].end = k;
            prev = k;
        }

        read.clear();

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        clock_gettime(CLOCK_MONOTONIC, &end);
        auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
        elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
        printf("(MSHLOADER)[LOAD]: End loading from parser, time elapsed: %f\n\n", elapsed);
#endif

        return 1;
    }


    static int dump(const Mesh<DIM, MESHSIZE> &mesh, VtkParser &parser,
                    const vector<double> &pointsData,
                    const vector<double> &elementData) {

        unordered_map<P, int> point_index;
        vector<array<double, 3>> points;
        int i = 0;

#ifdef MSHLOADER_VERBOSE
        timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        printf("\n(MSHLOADER)[DUMP]: Start dumping to parser\n");
        cout << "[DUMP]: starting dumping points ... ";
        cout.flush();
#endif

//
        for (const auto &point: mesh.points) {
            array<double, 3> temp = {0, 0, 0};
            for (int i = 0; i < DIM; ++i)
                temp[i] = point[i];
            points.push_back(temp);
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[DUMP]: starting dumping mesh elements ... ";
        cout.flush();
#endif
//
        vector<vector<int>> cells;
        for (const auto &tri_points: mesh.elements) {
            vector<int> t_cell;
            for (auto &p: tri_points)
                t_cell.emplace_back(p);
            cells.emplace_back(t_cell);
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[DUMP]: starting dumping points data ... ";
        cout.flush();
#endif

        vector<vector<double>> points_value;
        for (const auto &vals: pointsData) {
            vector<double> temp;
            temp.emplace_back(vals);
            points_value.emplace_back(temp);
        }
//        for (const auto &pair: mesh.index) {
//            if (pointsData.contains(pair.first)) {
//                vector<double> c_temp;
//                c_temp.emplace_back(pointsData.at(pair.first));
//                points_value.emplace_back(c_temp);
//            }
//        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
#endif

        vector<vector<double>> cell_value;

        parser.loadMesh(points, cells, points_value, cell_value);

#ifdef MSHLOADER_VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &end);
        auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
        elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
        printf("(MSHLOADER)[DUMP]: End dumping to parser, time elapsed: %f\n\n", elapsed);
#endif

        return 1;
    }
};


#endif //EIKONEL_TEST_MESHLOADER_H

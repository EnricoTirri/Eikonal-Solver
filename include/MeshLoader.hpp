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
struct MeshLoader {
    typedef MeshElement<DIM,MESHSIZE> M;
    typedef Eikonal_traits<DIM>::Point P;

    static int load(Mesh<DIM, MESHSIZE> &mesh, const VtkParser &parser,
                    unordered_map<P, double> &pointsData, unordered_map<M, double> &elementData) {
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
        unordered_map<P, vector<M *>> read;

#ifdef MSHLOADER_VERBOSE
        cout << "[LOAD]: starting loading points ... ";
        cout.flush();
#endif

        for (auto cell: parser.cells) {
            M &m = mesh.elements.emplace_back();
            for (int i = 0; i < MESHSIZE; ++i) {
                VtkPoint v_p = parser.points[cell.point_ids[i]];
                for (int j = 0; j < DIM; ++j)
                    m[i][j] = v_p.vec[j];
            }
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[LOAD]: starting loading mesh elements ... ";
        cout.flush();
#endif

        for (auto &element: mesh.elements) {
            for (int i = 0; i < MESHSIZE; i++) {
                read[element[i]].push_back(&element);
            }
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[LOAD]: starting loading adjacents ... ";
        cout.flush();
#endif

        size_t k = 0;
        size_t prev = 0;
        for (auto &i: read) {
            if (k == 0)
                mesh.index[i.first].start = 0;
            else
                mesh.index[i.first].start = prev;

            for (auto j: i.second) {
                mesh.adjacentList.push_back(j);
                k++;
            }
            mesh.index[i.first].end = k;
            prev = k;
        }

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
                    const unordered_map<P, double> &pointsData,
                    const unordered_map<M, double> &elementData) {

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

        for (const auto &pair: mesh.index) {
            point_index[pair.first] = i;
            array<double, 3> temp = {0, 0, 0};
            for (int i = 0; i < pair.first.size(); i++)
                temp[i] = pair.first[i];
            points.push_back(temp);
            ++i;
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[DUMP]: starting dumping mesh elements ... ";
        cout.flush();
#endif

        vector<vector<int>> cells;
        for (const auto &tri_points: mesh.elements) {
            vector<int> t_cell;
            for (auto &p: tri_points)
                t_cell.emplace_back(point_index[p]);
            cells.emplace_back(t_cell);
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[DUMP]: starting dumping points data ... ";
        cout.flush();
#endif

        vector<vector<double>> points_value;
        for (const auto &pair: mesh.index) {
            if (pointsData.contains(pair.first)) {
                vector<double> c_temp;
                c_temp.emplace_back(pointsData.at(pair.first));
                points_value.emplace_back(c_temp);
            }
        }

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

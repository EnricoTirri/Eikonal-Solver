
#include <MeshLoader.hpp>
#include <Mesh.h>
#include <vector>
#include <EikonalTraits.hpp>
#include <VtkParser.hpp>
#include <array>

#ifdef MSHLOADER_VERBOSE
#include <iostream>
#include <ctime>
using std::cout;
using std::endl;
#endif
namespace Eikonal {
    template<size_t MESHSIZE>
    int
    MeshLoader<MESHSIZE>::load(Mesh<MESHSIZE> &mesh, const VtkParser &parser, std::vector<double> &pointsData,
                                    std::vector<double> &elementData) {
        if (parser.status != 1)
            return 0;

#ifdef MSHLOADER_VERBOSE
        timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        printf("\n(MSHLOADER)[LOAD]: Start loading from parser\n");
#endif

        mesh.adjPointPtr.clear();
        mesh.pointAdjacentElementList.clear();
        mesh.adjElementPtr.clear();
        mesh.elementAdjacentPointList.clear();

#ifdef MSHLOADER_VERBOSE
        cout << "[LOAD]: starting loading points ... ";
        cout.flush();
#endif


        for (auto point: parser.points) {
            Traits::Point &p = mesh.points.emplace_back();
            for (int i = 0; i < 3; ++i)
                p[i] = point.vec[i];
        }

        mesh.adjElementPtr.reserve(parser.cells.size()+1);
        mesh.elementAdjacentPointList.resize(parser.cells.size() * MESHSIZE);
        int elPtr=0;
        for (auto cell: parser.cells) {
            mesh.adjElementPtr.push_back(elPtr);
#pragma unroll MESHSIZE
            for (int j = 0; j < MESHSIZE; ++j){
                mesh.elementAdjacentPointList[elPtr] = cell.point_ids[j];
                elPtr++;
            }
        }
        mesh.adjElementPtr.push_back(elPtr);

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[LOAD]: starting loading mesh elements_legacy ... ";
        cout.flush();
#endif


        //support structure for parsing file
        std::vector<std::vector<int>> read;
        read.resize(parser.points.size());

        for (int elID = 0; elID < mesh.adjElementPtr.size()-1; ++elID) {
            int sP = mesh.adjElementPtr[elID];
            for (int i = 0; i < MESHSIZE; i++) {
                int ptId = mesh.elementAdjacentPointList[sP+i];
                read[ptId].push_back(elID);
            }
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[LOAD]: starting loading adjacents ... ";
        cout.flush();
#endif

        size_t k = 0;
        mesh.adjPointPtr.resize(read.size() + 1);

        for (int i = 0; i < read.size(); ++i) {
            mesh.adjPointPtr[i] = k;

            for (auto j: read[i]) {
                mesh.pointAdjacentElementList.push_back(j);
                k++;
            }
        }

        mesh.adjPointPtr[read.size()] = k;

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

    template<size_t MESHSIZE>
    int MeshLoader<MESHSIZE>::dump(const Mesh<MESHSIZE> &mesh, VtkParser &parser,
                                        const std::vector<double> &pointsData,
                                        const std::vector<double> &elementData) {

        std::unordered_map<Traits::Point, int> point_index;
        std::vector<std::array<double, 3>> points;

#ifdef MSHLOADER_VERBOSE
        timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        printf("\n(MSHLOADER)[DUMP]: Start dumping to parser\n");
        cout << "[DUMP]: starting dumping points ... ";
        cout.flush();
#endif

//
        for (const auto &point: mesh.points) {
            std::array<double, 3> temp = {0, 0, 0};
            for (int i = 0; i < 3; ++i)
                temp[i] = point[i];
            points.push_back(temp);
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[DUMP]: starting dumping mesh elements_legacy ... ";
        cout.flush();
#endif
//
        std::vector<std::vector<int>> cells;
        for (int i=0; i<mesh.adjElementPtr.size()-1; ++i) {
            std::vector<int> t_cell;
            int start = mesh.adjElementPtr[i];
            int end = mesh.adjElementPtr[i+1];
            for (int j=start; j<end; ++j)
                t_cell.emplace_back(mesh.elementAdjacentPointList[j]);
            cells.emplace_back(t_cell);
        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
        cout << "[DUMP]: starting dumping points data ... ";
        cout.flush();
#endif

        std::vector<std::vector<double>> points_value;
        for (const auto &vals: pointsData) {
            std::vector<double> temp;
            temp.emplace_back(vals);
            points_value.emplace_back(temp);
        }
//        for (const auto &pair: mesh.adjPointPtr) {
//            if (pointsData.contains(pair.first)) {
//                std::vector<double> c_temp;
//                c_temp.emplace_back(pointsData.at(pair.first));
//                points_value.emplace_back(c_temp);
//            }
//        }

#ifdef MSHLOADER_VERBOSE
        cout << "end" << endl;
#endif

        std::vector<std::vector<double>> cell_value;

        parser.loadMesh(points, cells, points_value, cell_value);

#ifdef MSHLOADER_VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &end);
        auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
        elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
        printf("(MSHLOADER)[DUMP]: End dumping to parser, time elapsed: %f\n\n", elapsed);
#endif

        return 1;
    }
}

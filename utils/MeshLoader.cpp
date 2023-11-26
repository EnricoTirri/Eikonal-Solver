//
// Created by Enrico on 25/11/2023.
//

#include "MeshLoader.hpp"

int MeshLoader::load(Mesh &mesh, const VtkParser &parser, std::unordered_map<Point,double>& pointsData,
                     std::unordered_map<MeshElement,double>& elementData) {
    if (parser.status != 1)
        return 0;

    mesh.index.clear();
    mesh.elements.clear();
    mesh.adjacentList.clear();

    //support structure for parsing file
    std::unordered_map<Point, std::vector<MeshElement*>> read;

    for (auto cell: parser.cells) {
        MeshElement &m = mesh.elements.emplace_back();
        for(int i = 0; i<MESHSIZE; ++i){
            VtkPoint v_p = parser.points[cell.point_ids[i]];
            m[i] = Point{v_p.x, v_p.y};
        }
    }

    for(auto & element : mesh.elements){
        for(int i = 0; i<MESHSIZE; i++) {
            read[element[i]].push_back(&element);
        }
    }

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

    return 1;
}

int MeshLoader::dump(const Mesh &mesh, VtkParser &parser, const std::unordered_map<Point,double>& pointsData,
                     const std::unordered_map<MeshElement,double>& elementData) {

    std::unordered_map<Point, int> point_index;
    std::vector<std::array<double, 3>> points;
    int i = 0;
    for (const auto &pair: mesh.index) {
        point_index[pair.first] = i;
        points.push_back({pair.first[0], pair.first[1], 0});
        ++i;
    }

    std::vector<std::vector<int>> cells;
    for (const auto &tri_points: mesh.elements) {
        std::vector<int> t_cell;
        for (auto &p: tri_points)
            t_cell.emplace_back(point_index[p]);
        cells.emplace_back(t_cell);
    }

    std::vector<std::vector<double>> points_value;
    for (const auto& pair: mesh.index) {
        if(pointsData.contains(pair.first)) {
            std::vector<double> c_temp;
            c_temp.emplace_back(pointsData.at(pair.first));
            points_value.emplace_back(c_temp);
        }
    }

    std::vector<std::vector<double>> cell_value;

    parser.loadMesh(points, cells, points_value, cell_value);

    return 1;
}

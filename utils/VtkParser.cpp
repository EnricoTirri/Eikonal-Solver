//
// Created by Enrico on 23/11/2023.
//
#include <fstream>
#include <vector>
#include <array>
#include "VtkParser.hpp"
#include <iostream>

#ifndef PARSER_VERBOSE
#define PARSER_VERBOSE false
#endif

using namespace std;

void VtkParser::open(const std::string &filename) {
    ifstream in(filename);
    if(in.bad())
        return;

    getline(in, header);
    getline(in, description);
    getline(in, filetype);

    if (filetype != "ASCII") {
        std::cout << "error wrong file type" << std::endl;
        return;
    }


    ascii_parser(&in);
}

void VtkParser::ascii_parser(std::ifstream *in) {
    string temp;
    *in >> temp;
    if(temp!="DATASET")
        return;

    *in >> dataset_type;
    if(dataset_type!="UNSTRUCTURED_GRID") {
        printf("unsupported format");
        return;
    }

    bool is2DPoint = true;

    *in >> temp;
    while(!(*in).eof()){
        if (temp == "POINTS") {
            int n_points;
            *in >> n_points;

            *in >> temp; //reads datatype, discard, always use float

            points.clear();
#if PARSER_VERBOSE
            int step_point = n_points/10;
            int count = 0;
            cout << "(PARSER): starting parsing " << n_points << " points: ";
#endif
            for (int i = 0; i < n_points; ++i) {
#if PARSER_VERBOSE
                if(step_point == count){
                    cout << "0";
                    count = 0;
                }
                count++;
#endif
                double x, y, z;
                *in >> x;
                *in >> y;
                *in >> z;
                if(z>0) is2DPoint = false;
                points.emplace_back(x, y, z);
            }
            cout << endl;
        }else if(temp == "CELLS"){
            int n_cells;
            *in >> n_cells;

            *in >> temp; //read cell size, useless

#if PARSER_VERBOSE
            int step_cells = n_cells/10;
            int count = 0;
            cout << "(PARSER): starting parsing " << n_cells << " cells: ";
#endif
            for(int i = 0; i< n_cells; ++i){

#if PARSER_VERBOSE
                if(step_cells == count){
                    cout << "0";
                    count = 0;
                }
                count++;
#endif

                int n_points;
                *in >> n_points;

                VtkCell cell;
                for(int j = 0; j<n_points; ++j){
                    int id;
                    *in >> id;
                    cell.point_ids.emplace_back(id);
                }

                if(n_points>cell_max_d) cell_max_d = n_points;
                cells.emplace_back(cell);
            }
            cout << endl;
        }else if(temp == "CELL_TYPES"){
            int n_cells;
            *in >> n_cells;

#if PARSER_VERBOSE
            int step_cells = n_cells/10;
            int count = 0;
            cout << "(PARSER): starting parsing " << n_cells << " cells types: ";
#endif
            for(int i = 0; i<n_cells; ++i){
#if PARSER_VERBOSE
                if(step_cells == count){
                    cout << "0";
                    count = 0;
                }
                count++;
#endif

                int type;
                *in >> type;

                cells[i].type = type;
            }

#if PARSER_VERBOSE
            cout << endl;
#endif

        }
        *in >> temp;
    }

    is2DPoint ? point_max_d = 2 : point_max_d = 3;

    status = 1;
}


void VtkParser::save(const string &filename) {
    if(status != 1){
        printf("cannot save, parser not initialized correctly");
        return;
    }

    ofstream out(filename);

    out << header << endl;
    out << description << endl;
    out << filetype << endl;
    out << "DATASET " << dataset_type << endl;
    out << endl;

    out << "POINTS " << points.size() << " float" <<endl;
#if PARSER_VERBOSE
    cout << "(PARSER): starting saving points ... ";
#endif
    for(const auto& point : points){
        out << point.x << " " << point.y << " "<< point.z << endl;
    }
    out << endl;
#if PARSER_VERBOSE
    cout << "end" << endl;
    cout << "(PARSER): starting saving cells ... ";
#endif
    int cell_number_size = 0;
    for(const auto& cell : cells){
        ++cell_number_size;
        cell_number_size += cell.point_ids.size();
    }
    out << "CELLS " << cells.size() << " " << cell_number_size << endl;
    for(const auto& cell : cells){
        out << cell.point_ids.size();
        for(auto id : cell.point_ids){
            out << " " << id;
        }
        out << endl;
    }
    out << endl;

#if PARSER_VERBOSE
    cout << "end" << endl;
    cout << "(PARSER): starting saving cells types ... ";
#endif

    out << "CELL_TYPES " << cells.size() << endl;
    for(const auto& cell : cells){
        out << cell.type << endl;
    }
    out << endl;

#if PARSER_VERBOSE
    cout << endl;
#endif

    int point_data_size = 0;
    for(const auto& point: points){
        if(!point.data.empty()) ++point_data_size;
    }
    if(point_data_size!=0) {

#if PARSER_VERBOSE
        cout << "(PARSER): starting saving points data ... ";
#endif
        out << "POINT_DATA " << point_data_size << endl;
        out << "SCALARS points_table float 1" << endl;
        out << "LOOKUP_TABLE points_table" << endl;
        for (const auto &point: points) {
            for (auto value: point.data)
                out << " " << value;
            out << endl;
        }
        out << endl;

#if PARSER_VERBOSE
        cout << "end" << endl;
#endif
    }

    int cell_data_size = 0;
    for(const auto& cell: cells){
        if(!cell.data.empty()) ++cell_data_size;
    }
    if(cell_data_size>0) {
#if PARSER_VERBOSE
        cout << "(PARSER): starting saving cell data ... ";
#endif

        out << "CELL_DATA " << cell_data_size<< endl;
        out << "SCALARS cells_table float 1" << endl;
        out << "LOOKUP_TABLE cells_table" << endl;
        for (const auto &cell: cells) {
            for (auto value: cell.data)
                out << " " << value;
            out << endl;
        }

#if PARSER_VERBOSE
        cout << "end" << endl;
#endif
    }
}

void VtkParser::loadMesh(const std::vector<std::array<double,3>>& new_points, const std::vector<std::vector<int>>& new_cells,
                     const std::vector<std::vector<double>>& point_data, const std::vector<std::vector<double>>& cell_data) {
    status = 0;
    points.clear();
    cells.clear();

#if PARSER_VERBOSE
    cout << "(PARSER): starting loading points and points data ... ";
#endif
    for(int i = 0; i<new_points.size();++i){
        std::vector<double> data;
        if(i<point_data.size())
            data = point_data[i];
        points.emplace_back(new_points[i][0], new_points[i][1], new_points[i][2], data);
    }

#if PARSER_VERBOSE
    cout << "end" << endl;
    cout << "(PARSER): starting loading cells and cell data ... ";
#endif

    for(int i = 0; i<new_cells.size();++i){
        std::vector<double> data;
        if(i<cell_data.size())
            data = cell_data[i];
        cells.emplace_back(5,new_cells[i],data);
    }

#if PARSER_VERBOSE
    cout << "end" << endl;
#endif

    dataset_type = "UNSTRUCTURED_GRID";
    header = "# vtk DataFile Version 2.0";
    description = "triangular elements loaded from computed data";
    filetype = "ASCII";
    status = 1;
}


//
// Created by Enrico on 23/11/2023.
//
#include <fstream>
#include <vector>
#include <hash_map>
#include "VtkParser.hpp"

using namespace std;

void VtkParser::open(const std::string &filename) {
    ifstream in(filename);
    if(in.bad())
        return;

    getline(in, header);
    getline(in, description);
    getline(in, filetype);

    if(filetype!="ASCII")
        return;

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

    *in >> temp;
    while(!(*in).eof()){
        if (temp == "POINTS") {
            int n_points;
            *in >> n_points;

            *in >> temp; //reads datatype, discard, always use float

            points.clear();
            for (int i = 0; i < n_points; ++i) {
                float x, y, z;
                *in >> x;
                *in >> y;
                *in >> z;
                points.emplace_back(x, y, z);
            }
        }else if(temp == "CELLS"){
            int n_cells;
            *in >> n_cells;

            *in >> temp; //read cell size, useless

            for(int i = 0; i< n_cells; ++i){
                int n_points;
                *in >> n_points;

                VtkCell cell;
                for(int j = 0; j<n_points; ++j){
                    int id;
                    *in >> id;
                    cell.point_ids.emplace_back(id);
                }

                cells.emplace_back(cell);
            }
        }else if(temp == "CELL_TYPES"){
            int n_cells;
            *in >> n_cells;

            for(int i = 0; i<n_cells; ++i){
                u_int8_t type;
                *in >> type;

                cells[i].type = type;
            }
        }
        *in >> temp;
    }
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
    for(auto point : points){
        out << point.x << " " << point.y << " "<< point.z << endl;
    }
    out << endl;

    int datasize = 0;
    for(auto cell : cells){
        ++datasize;
        datasize += cell.point_ids.size();
    }
    out << "CELLS " << cells.size() << " " << datasize << endl;
    for(auto cell : cells){
        out << cell.point_ids.size();
        for(auto id : cell.point_ids){
            out << " " << id;
        }
        out << endl;
    }
    out << endl;

    out << "CELL_TYPES " << cells.size() << endl;
    for(auto cell : cells){
        out << cell.type << endl;
    }
}


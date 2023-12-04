//
// Created by Enrico on 23/11/2023.
//
#include <fstream>
#include <utility>

#ifndef EIKONEL_TEST_MSHPARSER_H
#define EIKONEL_TEST_MSHPARSER_H

class VtkPoint{
public:
    std::array<double,3> vec;
    std::vector<double> data;

    VtkPoint(double x, double y, double z):vec{x,y,z}{}
    VtkPoint(double x, double y, double z, std::vector<double> data) : VtkPoint(x,y,z){
        this->data = std::move(data);
    };


    double x() const { return vec[0]; }

    double y() const { return vec[1]; }

    double z() const { return vec[2]; }
};


class VtkCell{
public:
    int type;
    std::vector<double> data;
    std::vector<int> point_ids;
    VtkCell() : type(-1){};
    VtkCell(int type, std::vector<int> point_ids) : type(type) { this->point_ids = std::move(point_ids); }
    VtkCell(int type, std::vector<int> point_ids, std::vector<double> data) : VtkCell(type,std::move(point_ids)){
        this->data = std::move(data);
    }


};


class VtkParser {
public:
    //status of parsing: 1=successful, 0=not successful
    int status = 0;

    //elements format info
    std::string header;
    std::string description;
    std::string filetype;

    std::string dataset_type;

    //nodes info
    std::vector<VtkPoint> points;

    //elements info
    std::vector<VtkCell> cells;

    void loadMesh(const std::vector<std::array<double,3>>& points,
              const std::vector<std::vector<int>>& cells,
              const std::vector<std::vector<double>>& point_data,
              const std::vector<std::vector<double>>& cell_data);

    void open(std::string const &filename);

    void save(std::string const &filename);


private:
    void ascii_parser(std::ifstream  &input);
    void ascii_saver(std::ofstream  &output);
};

#endif

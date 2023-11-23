//
// Created by Enrico on 23/11/2023.
//
#include <fstream>

#ifndef EIKONEL_TEST_MSHPARSER_H
#define EIKONEL_TEST_MSHPARSER_H

class VtkPoint{
public:
    float x,y,z;
    VtkPoint(float x, float y, float z):x{x},y{y},z{z}{}
};


class VtkCell{
public:
    u_int8_t type;
    std::vector<int> point_ids;
};


class VtkParser {
public:
    //status of parsing: 1=successful, 0=not successful
    int status;

    //mesh format info
    std::string header;
    std::string description;
    std::string filetype;

    std::string dataset_type;

    //nodes info
    std::vector<VtkPoint> points;

    //elements info
    std::vector<VtkCell> cells;


    VtkParser(){
        status = 0;
        header = "";
        description = "";
        filetype = "";
        dataset_type = "";
    }

    void open(std::string const &filename);

    void save(std::string const &filename);

private:
    void ascii_parser(std::ifstream *input);
};
#endif //EIKONEL_TEST_MSHPARSER_H

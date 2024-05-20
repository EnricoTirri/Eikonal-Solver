
#ifndef VTK_PARSER
#define VTK_PARSER

#include <fstream>
#include <utility>
#include <array>
#include <vector>

class VtkPoint {
public:
    std::array<double, 3> vec;
    std::vector<double> data;

    VtkPoint(double x, double y, double z) : vec{x, y, z} {}

    VtkPoint(double x, double y, double z, std::vector<double> data) : VtkPoint(x, y, z) {
        this->data = std::move(data);
    };


    double x() const { return vec[0]; }

    double y() const { return vec[1]; }

    double z() const { return vec[2]; }
};


class VtkCell {
public:
    int type;
    std::vector<double> data;
    std::vector<int> point_ids;

    VtkCell() : type(-1) {};

    VtkCell(int type, std::vector<int> point_ids) : type(type) { this->point_ids = std::move(point_ids); }

    VtkCell(int type, std::vector<int> point_ids, std::vector<double> data) : VtkCell(type, std::move(point_ids)) {
        this->data = std::move(data);
    }
};


class VtkParser {
public:
    //status of parsing: 1=successful, 0=not successful
    int status = 0;

    //elements_legacy format info
    std::string header;
    std::string description;
    std::string filetype;

    std::string dataset_type;

    //nodes info
    std::vector<VtkPoint> points;

    //elements_legacy info
    std::vector<VtkCell> cells;

    void loadMesh(const std::vector<std::array<double, 3>> &points,
                  const std::vector<std::vector<int>> &cells,
                  const std::vector<std::vector<double>> &point_data,
                  const std::vector<std::vector<double>> &cell_data);

    void open(std::string const &filename);

    void save(std::string const &filename);


private:
    void ascii_parser(std::ifstream &input);

    void ascii_saver(std::ofstream &output);
};

#endif //VTK_PARSER
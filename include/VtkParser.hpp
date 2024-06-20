
#ifndef VTK_PARSER
#define VTK_PARSER

#include <fstream>
#include <utility>
#include <array>
#include <vector>

// This class describe the data structure of a point parsed from .vtk file
class VtkPoint {
public:
    // coordinates of the point
    std::array<double, 3> vec;
    // data associated to the point
    std::vector<double> data;

    VtkPoint(double x, double y, double z) : vec{x, y, z} {}

    VtkPoint(double x, double y, double z, std::vector<double> data) : VtkPoint(x, y, z) {
        this->data = std::move(data);
    };


    double x() const { return vec[0]; }

    double y() const { return vec[1]; }

    double z() const { return vec[2]; }
};

// This class describe the data structure of a cell parsed from .vtk file
class VtkCell {
public:
    // Type of cell
    int type;
    // data associated to the cell
    std::vector<double> data;
    // ids of points the cell is composed by
    std::vector<int> point_ids;

    VtkCell() : type(-1) {};

    VtkCell(int type, std::vector<int> point_ids) : type(type) { this->point_ids = std::move(point_ids); }

    VtkCell(int type, std::vector<int> point_ids, std::vector<double> data) : VtkCell(type, std::move(point_ids)) {
        this->data = std::move(data);
    }
};

// This class represent the data structure of a .vtk file
// contains method to read and write .vtk files
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

    // This method initialize the class using the provided data
    void loadMesh(const std::vector<std::array<double, 3>> &points,
                  const std::vector<std::vector<int>> &cells,
                  const std::vector<std::vector<double>> &point_data,
                  const std::vector<std::vector<double>> &cell_data);

    // This method reads a .vtk file loading into class data
    void open(std::string const &filename);

    // This method write class data to a .vtk file
    void save(std::string const &filename);


private:
    // This method contains algorithm to read a ASCII formatted .vtk file
    void ascii_parser(std::ifstream &input);

    // This method contains algorithm to write a ASCII formatted .vtk file
    void ascii_saver(std::ofstream &output);
};

#endif //VTK_PARSER
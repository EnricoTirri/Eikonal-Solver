//
// Created by Enrico on 23/11/2023.
//
#include <fstream>
#include <vector>
#include <array>
#include "VtkParser.hpp"
#include <iostream>

using namespace std;

void VtkParser::open(const std::string &filename) {
    ifstream in(filename);
    if (in.bad())
        return;

    getline(in, header);
    getline(in, description);
    getline(in, filetype);


#ifdef PARSER_VERBOSE
    timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    printf("\n(PARSER)[OPEN]: Start parsing from file\n");
#endif

    ascii_parser(in);

#ifdef PARSER_VERBOSE
    clock_gettime(CLOCK_MONOTONIC, &end);
    auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
    elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
    printf("(PARSER)[OPEN]: End parsing from file, time elapsed: %f\n\n", elapsed);
#endif
}

void VtkParser::ascii_parser(std::ifstream &in) {
    string temp;
    in >> temp;
    if (temp != "DATASET")
        return;

    in >> dataset_type;
    if (dataset_type != "UNSTRUCTURED_GRID") {
        printf("[OPEN]: unsupported format");
        return;
    }

    in >> temp;
    while (!(in).eof()) {
        if (temp == "POINTS") {
            int n_points;
            in >> n_points;

            in >> temp; //reads datatype, discard, always use float

            points.clear();
#ifdef PARSER_VERBOSE
            int step_point = (n_points / 10) - 1;
            int step = 0;
            int count = 0;
            cout << "[OPEN]: starting parsing " << n_points << " points: ";
#endif
            for (int i = 0; i < n_points; ++i) {
#ifdef PARSER_VERBOSE
                if (step_point == count) {
                    step += 10;
                    cout << step << "% ";
                    cout.flush();
                    count = 0;
                }
                count++;
#endif
                double x, y, z;
                in >> x;
                in >> y;
                in >> z;
                points.emplace_back(x, y, z);
            }
#ifdef PARSER_VERBOSE
            cout << endl;
#endif
        } else if (temp == "CELLS") {
            int n_cells;
            in >> n_cells;

            in >> temp; //read cell size, useless

#ifdef PARSER_VERBOSE
            int step_cells = (n_cells / 10) - 1;
            int count = 0;
            int step = 0;
            cout << "[OPEN]: starting parsing " << n_cells << " cells: ";
#endif
            for (int i = 0; i < n_cells; ++i) {

#ifdef PARSER_VERBOSE
                if (step_cells == count) {
                    step += 10;
                    cout << step << "% ";
                    cout.flush();
                    count = 0;
                }
                count++;
#endif

                int n_points;
                in >> n_points;

                VtkCell cell;
                for (int j = 0; j < n_points; ++j) {
                    int id;
                    in >> id;
                    cell.point_ids.emplace_back(id);
                }
                cells.emplace_back(cell);
            }
#ifdef PARSER_VERBOSE
            cout << endl;
#endif
        } else if (temp == "CELL_TYPES") {
            int n_cells;
            in >> n_cells;

#ifdef PARSER_VERBOSE
            int step_cells = (n_cells / 10) - 1;
            int step = 0;
            int count = 0;
            cout << "[OPEN]: starting parsing " << n_cells << " cells types: ";
            cout.flush();
#endif
            for (int i = 0; i < n_cells; ++i) {
#ifdef PARSER_VERBOSE
                if (step_cells == count) {
                    step += 10;
                    cout << step << "% ";
                    cout.flush();
                    count = 0;
                }
                count++;
#endif

                int type;
                in >> type;

                cells[i].type = type;
            }

#ifdef PARSER_VERBOSE
            cout << endl;
#endif

        }
        in >> temp;
    }

    status = 1;
}

void VtkParser::save(const std::string &filename) {
    ofstream out(filename);
    if (out.bad())
        return;

#ifdef PARSER_VERBOSE
    timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    printf("\n(PARSER)[SAVE]: Start saving to file\n");
#endif

    ascii_saver(out);

#ifdef PARSER_VERBOSE
    clock_gettime(CLOCK_MONOTONIC, &end);
    auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
    elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
    printf("(PARSER)[SAVE]: End saving to file, time elapsed: %f\n\n", elapsed);
#endif
}

void VtkParser::ascii_saver(ofstream &out) {
    if (status != 1) {
        printf("[SAVE]: Cannot save, parser not initialized correctly\n");
        return;
    }

    out << header << endl;
    out << description << endl;
    out << filetype << endl;
    out << "DATASET " << dataset_type << endl;
    out << endl;

    out << "POINTS " << points.size() << " float" << endl;
#ifdef PARSER_VERBOSE
    cout << "[SAVE]: starting saving points ... ";
    cout.flush();
#endif
    for (auto &point: points) {
        out << point.x() << " " << point.y() << " " << point.z() << endl;
    }
    out << endl;
#ifdef PARSER_VERBOSE
    cout << "end" << endl;
    cout << "[SAVE]: starting saving cells ... ";
    cout.flush();
#endif
    int cell_number_size = 0;
    for (const auto &cell: cells) {
        ++cell_number_size;
        cell_number_size += cell.point_ids.size();
    }
    out << "CELLS " << cells.size() << " " << cell_number_size << endl;
    for (const auto &cell: cells) {
        out << cell.point_ids.size();
        for (auto id: cell.point_ids) {
            out << " " << id;
        }
        out << endl;
    }
    out << endl;

#ifdef PARSER_VERBOSE
    cout << "end" << endl;
    cout << "[SAVE]: starting saving cells types ... ";
    cout.flush();
#endif

    out << "CELL_TYPES " << cells.size() << endl;
    for (const auto &cell: cells) {
        out << cell.type << endl;
    }
    out << endl;

#ifdef PARSER_VERBOSE
    cout << "end" << endl;
#endif

    int point_data_size = 0;
    for (const auto &point: points) {
        if (!point.data.empty()) ++point_data_size;
    }
    if (point_data_size != 0) {

#ifdef PARSER_VERBOSE
        cout << "[SAVE]: starting saving points data ... ";
        cout.flush();
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

#ifdef PARSER_VERBOSE
        cout << "end" << endl;
#endif
    }

    int cell_data_size = 0;
    for (const auto &cell: cells) {
        if (!cell.data.empty()) ++cell_data_size;
    }
    if (cell_data_size > 0) {

#ifdef PARSER_VERBOSE
        cout << "[SAVE]: starting saving cell data ... ";
        cout.flush();
#endif

        out << "CELL_DATA " << cell_data_size << endl;
        out << "SCALARS cells_table float 1" << endl;
        out << "LOOKUP_TABLE cells_table" << endl;
        for (const auto &cell: cells) {
            for (auto value: cell.data)
                out << " " << value;
            out << endl;
        }

#ifdef PARSER_VERBOSE
        cout << "end" << endl;
#endif
    }
}

void VtkParser::loadMesh(const std::vector<std::array<double, 3>> &new_points,
                         const std::vector<std::vector<int>> &new_cells,
                         const std::vector<std::vector<double>> &point_data,
                         const std::vector<std::vector<double>> &cell_data) {
    status = 0;
    points.clear();
    cells.clear();

#ifdef PARSER_VERBOSE
    timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    printf("\n(PARSER)[LOAD]: Start loading from data\n");
    cout << "[LOAD]: starting loading points and points data ... ";
    cout.flush();
#endif
    for (int i = 0; i < new_points.size(); ++i) {
        std::vector<double> data;
        if (i < point_data.size())
            data = point_data[i];
        points.emplace_back(new_points[i][0], new_points[i][1], new_points[i][2], data);
    }

#ifdef PARSER_VERBOSE
    cout << "end" << endl;
    cout << "[LOAD]: starting loading cells and cell data ... ";
    cout.flush();
#endif

    for (int i = 0; i < new_cells.size(); ++i) {
        std::vector<double> data;
        if (i < cell_data.size())
            data = cell_data[i];
        cells.emplace_back(5, new_cells[i], data);
    }

#ifdef PARSER_VERBOSE
    cout << "end" << endl;
    clock_gettime(CLOCK_MONOTONIC, &end);
    auto elapsed = static_cast<double>((end.tv_sec - start.tv_sec));
    elapsed += static_cast<double>((end.tv_nsec - start.tv_nsec)) / 1000000000.0;
    printf("(PARSER)[LOAD]: End loading from data, time elapsed: %f\n\n", elapsed);
#endif

    dataset_type = "UNSTRUCTURED_GRID";
    header = "# vtk DataFile Version 2.0";
    description = "triangular elements loaded from computed data";
    filetype = "ASCII";
    status = 1;
}


#include "LocalProblem/include/SimplexData.hpp"
#include "LocalProblem/include/solveEikonalLocalProblem.hpp"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include "VtkParser.hpp"

#ifndef TYPE
#define TYPE 2
#endif


#include <vector>
#include <Eigen/Core>
#include <PointsEdge.h>

#define MAXFLOAT 900000
//now we will implement this algorithm Fast iterative method (X,L)
//define hash function for Point


/*double solveQuadratic(double a, double b, double c){
    double u = c + 1;
    if (u <= b) return u;
    u = (b+c+std::sqrt(-(b*b) - (c*c) + 2*b*c +2))/2;
    if(u <= a) return u;
    u = (2*(a+b+c) + std::sqrt(4*(a+b+c)*(a+b+c) - 12*(a*a + b*b + c*c -1)))/6;
    return u;
}*/

//X will be the starting point from witch the wave will propagate, L is the active list
void FIM(std::unordered_map<Point, double> &U, std::vector<Point> X, std::vector<Point> L, PointsEdge &data, FILE *fp) {

    U.reserve(data.index.size());

//#pragma omp parallel for
    for (auto &i: data.index) {
        U.insert({i.first, MAXFLOAT});
    }
    for (const auto &i: X) {
        U[i] = 0;
        L.push_back(i);
    }

    //2. Update points in L
    while (!L.empty()) {
        std::vector<Point> next_L;
        next_L.clear();
//#pragma omp parallel for//matbe task would be better
        for (const auto &i: L) {
            //take time value of point p
            double p = U[i];
            //find neighbors of L[i] and get the base (the DIMENSION points with the smallest value of U
            std::vector<mesh_element> neighbors;
            std::size_t start, end;
            start = data.index[i].start;
            end = data.index[i].end;
            for (std::size_t j = start; j < end; j++) {
                neighbors.push_back(static_cast<mesh_element>(*data.adjacentList[j]));
            }

            //with task maybe we can parallelize this
            for (const auto &m_element: neighbors) {
                for (auto &point: m_element) {
                    if (point == i) continue;
                    //if point in L continue
                    //solve local problem with this point as unknown and the others as base
                    std::array<Point, DIMENSION + 1> base;
                    if (point[0] == 2. && point[1] == 2) {
                        for (int j = 0; j < DIMENSION + 1; ++j) {
                            printf("x: %f y: %f \t U: %f\n", m_element[j][0], m_element[j][1], U[m_element[j]]);
                        }

                        printf("found point\n");
                    }

                    std::size_t k = 0;
                    Eikonal::Eikonal_traits<DIMENSION>::VectorExt values;
                    for (std::size_t j = 0; j < DIMENSION + 1; j++) {
                        if (m_element[j] == point) {
                            base[DIMENSION] = point;
                        } else {
                            base[k] = m_element[j];
//                            values[static_cast<int>(k)] = U[m_element[j]];
                            k++;
                        }
                    }
                    for (int iter = 0; iter < DIMENSION; iter++) {
                        values[iter] = U[base[iter]];
                    }

                    Eikonal::Eikonal_traits<DIMENSION>::MMatrix M;
                    M << 1.0, 0.0,
                            0.0, 1.0;
//                    printf("base: ");
//                    for (int j = 0; j < DIMENSION +1; ++j) {
//                        printf("x: %f y: %f \t",base[j][0],base[j][1]);
//                    }
//                    printf("\n while examining point x: %f y:%f\n",point[0],point[1]);
                    Eikonal::SimplexData<DIMENSION> simplex{base, M};
                    Eikonal::solveEikonalLocalProblem<DIMENSION> solver{std::move(simplex),
                                                                        values};
                    auto sol = solver();
                    //if no descent direction or no convergence kill the process
                    if (sol.status != 0) {
                        printf("error on convergence\n");
                        return;
                    }
                    auto newU = sol.value;
/*
                    double b=-1;
                    double c=-1;
                    for(Point i2 : base)
                    {
                        if(i2!=point)
                        {
                            if (c<0)
                                c=U[i2];
                            else
                                b=U[i2];
                        }
                    }
                    auto newU2 = solveQuadratic(99,b,c);
                    if (abs(newU2-newU)>0.1)
                    {
                       // printf("difference local %f, quad: %f old: %f\n",newU,newU2,U[point]);
                    }
                    newU = newU2;*/
                    if (newU < U[point]) {
                        U[point] = newU;
                        // #pragma omp atomic
                        next_L.push_back(point);
                    }
                }
            }
        }
        //remove duplicates from next_L
        std::sort(next_L.begin(), next_L.end(), [](const Point &a, const Point &b) {
            return a[0] * 1000 + a[1] < b[0] * 1000 + b[1];
        });
        next_L.erase(std::unique(next_L.begin(), next_L.end()), next_L.end());
        for (auto &i: L) {
            fprintf(fp, "%f %f %f\n", i[0], i[1], U[i]);
        }
        fprintf(fp, "end queue\n");
        L = std::move(next_L);
    }
}


int main() {
    VtkParser parser;
    parser.open("testmesh.vtk");

    //support structure for parsing file
    std::unordered_map<Point, std::vector<mesh_element>> readed;

    //final structure
    PointsEdge pointsEdge;


    std::vector<mesh_element> temp_mesh_list;
    for(auto cell : parser.cells){
        Point a = Point {parser.points[cell.point_ids[0]].x,parser.points[cell.point_ids[0]].y};
        Point b = Point {parser.points[cell.point_ids[1]].x,parser.points[cell.point_ids[1]].y};
        Point c = Point {parser.points[cell.point_ids[2]].x,parser.points[cell.point_ids[2]].y};
        temp_mesh_list.emplace_back(mesh_element{a, b, c});
    }

    for(auto mesh : temp_mesh_list){
        readed[mesh[0]].emplace_back(mesh);
        readed[mesh[1]].emplace_back(mesh);
        readed[mesh[2]].emplace_back(mesh);
        pointsEdge.mesh.emplace_back(mesh);
    }

    size_t k = 0;
    size_t prev = 0;
    for (auto &i: readed) {
        if (k == 0)
            pointsEdge.index[i.first].start = 0;
        else
            pointsEdge.index[i.first].start = prev;

        for (auto &j: i.second) {
            pointsEdge.adjacentList.push_back(&j);
            k++;
        }
        pointsEdge.index[i.first].end = k;
        prev = k;
    }
    printf("end parsing\n");



    //let's try it
    std::unordered_map<Point, double> U;
    std::vector<Point> L;
    std::vector<Point> X;

    X.emplace_back(0, 3);
    X.emplace_back(10, 3);
    X.emplace_back(0, 6);
    X.emplace_back(10, 6);
    time_t start = clock();
    FILE *fp = fopen("outL.txt", "w");
    FIM(U, X, L, pointsEdge, fp);
    time_t end = clock();
    printf("time elapsed: %f\n", (double) (end - start) / CLOCKS_PER_SEC);
    printf("end FIM\n");


    std::unordered_map<Point,int> point_index;
    std::vector<std::array<double,3>> points;
    int i = 0;
    for(const auto& pair : pointsEdge.index){
        point_index[pair.first] = i;
        points.push_back({pair.first[0], pair.first[1], 0});
        ++i;
    }

    std::vector<std::vector<int>> cells;
    for(const auto& tri_points : pointsEdge.mesh){
        std::vector<int> t_cell;
        for(auto &p : tri_points)
            t_cell.emplace_back(point_index[p]);
        cells.emplace_back(t_cell);
    }

    std::vector<std::vector<double>> points_value;
    for(auto pair : pointsEdge.index){
        std::vector<double> c_temp;
        c_temp.emplace_back(U[pair.first]);
        points_value.emplace_back(c_temp);
    }

    std::vector<std::vector<double>> cell_value;

    parser.loadTriangular(points, cells, points_value, cell_value);

    parser.save("final_out.vtk");

    return 0;
}

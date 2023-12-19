//
// Created by Enrico on 19/12/2023.
//
#include <iostream>
#include <vector>
#include <Mesh.h>

#ifndef EIKONAL_SOLVER_EIKONALSOLVER_H
#define EIKONAL_SOLVER_EIKONALSOLVER_H

template<int DIM,int MESH_SIZE>
class EikonalSolver{
public:
    double const MAXF = 900000;

    bool solve(std::vector<double> &U,
                        const std::vector<int>& X,
                        const Mesh<DIM, MESH_SIZE> &data);

    void print_spec();
};

template class EikonalSolver<2,3>;
template class EikonalSolver<3,3>;
template class EikonalSolver<3,4>;

#endif //EIKONAL_SOLVER_EIKONALSOLVER_H

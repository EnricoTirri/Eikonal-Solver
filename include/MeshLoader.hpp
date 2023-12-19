//
// Created by Enrico on 25/11/2023.
//
#include <Mesh.h>
#include <Eikonal_traits.hpp>
#include <VtkParser.hpp>
#include <vector>

template<size_t DIM, size_t MESHSIZE>
class MeshLoader{

    typedef Mesh<DIM,MESHSIZE>::MeshElement M;
    typedef Eikonal_traits<DIM>::Point P;

public:
    int load(Mesh<DIM, MESHSIZE> &mesh, const VtkParser &parser,
                    std::vector<double> &pointsData, std::vector<double> &elementData);


    int dump(const Mesh<DIM, MESHSIZE> &mesh, VtkParser &parser,
                    const std::vector<double> &pointsData,
                    const std::vector<double> &elementData);
};

template class MeshLoader<2,3>;
template class MeshLoader<3,3>;
template class MeshLoader<3,4>;


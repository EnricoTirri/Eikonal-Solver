//
// Created by Enrico on 25/11/2023.
//

#ifndef EIKONEL_TEST_MESHLOADER_H
#define EIKONEL_TEST_MESHLOADER_H

#include <Mesh.h>
#include <Eikonal_traits.hpp>
#include <VtkParser.hpp>

using namespace Eikonal;
using namespace std;

template<size_t DIM, size_t MESHSIZE>
class MeshLoader {
    typedef MeshElement<DIM,MESHSIZE> M;
    typedef Eikonal_traits<DIM>::Point P;

public:
    static int load(Mesh<DIM, MESHSIZE> &mesh, const VtkParser &parser,
                    unordered_map<P, double> &pointsData,
                    unordered_map<M, double> &elementData);

    static int dump(const Mesh<DIM,MESHSIZE> &mesh, VtkParser &parser,
                    const unordered_map<P, double> &pointsData,
                    const unordered_map<M, double> &elementData);
};


#endif //EIKONEL_TEST_MESHLOADER_H

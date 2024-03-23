//
// Created by Enrico on 25/11/2023.
//
#include <Mesh.h>
#include <EikonalTraits.hpp>
#include <VtkParser.hpp>
#include <vector>

namespace Eikonal {
    template< size_t MESHSIZE>
    class MeshLoader {

        typedef Mesh<MESHSIZE>::MeshElement M;

    public:
        int load(Mesh<MESHSIZE> &mesh, const VtkParser &parser,
                 std::vector<double> &pointsData, std::vector<double> &elementData);


        int dump(const Mesh<MESHSIZE> &mesh, VtkParser &parser,
                 const std::vector<double> &pointsData,
                 const std::vector<double> &elementData);
    };

    template
    class MeshLoader<3>;

    template
    class MeshLoader<4>;

}
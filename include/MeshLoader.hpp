
#ifndef MESH_LOADER
#define MESH_LOADER

#include "Mesh.hpp"
#include "VtkParser.hpp"
#include <vector>

namespace Eikonal {
    template<size_t MESHSIZE>
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

#endif //MESH_LOADER

#ifndef MESH_LOADER
#define MESH_LOADER

#include "Mesh.hpp"
#include "VtkParser.hpp"
#include <vector>

namespace Eikonal {

    // This class defines the methods to transfer Mesh data from Parser to Mesh data structure
    template<size_t MESHSIZE>
    class MeshLoader {

        typedef Mesh<MESHSIZE>::MeshElement M;

    public:
        // This method load Mesh data to Mesh data structure, return 1 if loaded correctly
        int load(Mesh<MESHSIZE> &mesh, const VtkParser &parser,
                 std::vector<double> &pointsData, std::vector<double> &elementData);

        // This method dump data to Parser, return 1 if dumped correctly
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
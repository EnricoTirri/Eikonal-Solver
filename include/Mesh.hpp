
#ifndef MESH
#define MESH

#include <vector>
#include "EikonalTraits.hpp"

namespace Eikonal {

    // This class defines the data structure for the Mesh
    template<std::size_t MESHSIZE>
    class Mesh {
    public:
        // Mesh element typedef
        using MeshElement = std::array<int, MESHSIZE>;

        // List of points of Mesh
        std::vector<typename Traits::Point> points;

        // CSR data structure Point to Element
        std::vector<int32_t> pointAdjacentElementList;
        std::vector<int32_t> adjPointPtr;

        // CSR data structure Element to Points
        std::vector<int32_t> adjElementPtr;
        std::vector<int32_t> elementAdjacentPointList;
    };
}

#endif //MESH
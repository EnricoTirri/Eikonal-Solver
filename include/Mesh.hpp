
#ifndef MESH
#define MESH

#include <vector>
#include "EikonalTraits.hpp"

namespace Eikonal {
    template<std::size_t MESHSIZE>
    class Mesh {
    public:
        using MeshElement = std::array<int, MESHSIZE>;

        std::vector<int32_t> pointAdjacentElementList;
        std::vector<int32_t> adjPointPtr;
        std::vector<typename Traits::Point> points;
        std::vector<int32_t> adjElementPtr;
        std::vector<int32_t> elementAdjacentPointList;
    };
}

#endif //MESH
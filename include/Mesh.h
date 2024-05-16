//
// Created by giorgio on 20/11/23.
//

#ifndef EIKONEL_TEST_MESH_H
#define EIKONEL_TEST_MESH_H

#include <Eigen/Core>
#include <vector>
#include <memory>
#include <EikonalTraits.hpp>

namespace Eikonal {
    template<std::size_t MESHSIZE>
    class Mesh {
    public:
        using MeshElement = std::array<int, MESHSIZE>;

        std::vector<int32_t> adjacentList;
        std::vector<int32_t> index;
        std::vector<typename Traits::Point> points;
        std::vector<MeshElement> elements_legacy;
    };
}

#endif

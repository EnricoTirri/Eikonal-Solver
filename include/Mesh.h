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
    typedef struct {
        size_t start, end;
    } range_t;

    template<std::size_t MESHSIZE>
    class Mesh {
    public:
        using MeshElement = std::array<int, MESHSIZE>;

        std::vector<uint32_t> adjacentList;
        std::vector<range_t> index;
        std::vector<typename Traits::Point> points;
        std::vector<MeshElement> elements_legacy;
        std::vector<uint32_t> eind;
        std::vector<uint32_t> eptr;
    };
}

#endif

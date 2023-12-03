//
// Created by Enrico on 28/11/2023.
//

#include <Eikonal_traits.hpp>

namespace Eikonal {
    template<std::size_t DIM, std::size_t MESHSIZE>
    using MeshElement = std::array<int, MESHSIZE>;
}
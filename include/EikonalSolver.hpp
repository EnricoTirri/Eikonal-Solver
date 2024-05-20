
#ifndef EIKONAL_SOLVER
#define EIKONAL_SOLVER

#include <vector>
#include "Mesh.hpp"

namespace Eikonal {
    template<int MESH_SIZE>
    class EikonalSolver {
    public:
        constexpr static double const MAXF = 900000;

        bool solve(std::vector<double> &U,
                   const std::vector<int> &X,
                   const Mesh<MESH_SIZE> &data);

        void print_spec();
    };

    template
    class EikonalSolver<3>;

    template
    class EikonalSolver<4>;
}

#endif //EIKONAL_SOLVER
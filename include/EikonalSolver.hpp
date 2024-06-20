
#ifndef EIKONAL_SOLVER
#define EIKONAL_SOLVER

#include <vector>
#include "Mesh.hpp"

// This class defines the structure of a Eikonal Global problem solver
namespace Eikonal {
    template<int MESH_SIZE>
    class EikonalSolver {
    public:
        // Time references of compute_time and prepare_time
        double compute = -1;
        double prepare = -1;

        // Reference value for initialization of solution
        constexpr static double const MAXF = 900000;

        // This method performs the global solver algorithms on Mesh "data",
        // using starting points "X" and returning solution vector "U"
        bool solve(std::vector<double> &U,
                   const std::vector<int> &X,
                   const Mesh<MESH_SIZE> &data);

        // This method prints general specification info of the solver
        void print_spec();
    };

    template
    class EikonalSolver<3>;


    template
    class EikonalSolver<4>;
}

#endif //EIKONAL_SOLVER
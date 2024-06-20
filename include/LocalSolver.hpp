
#ifndef LOCAL_SOLVER
#define LOCAL_SOLVER

#include "EikonalTraits.hpp"

namespace Eikonal {

    // This class defines the structure of the Eikonal Local solver
    template<int MESH_SIZE>
    class LocalSolver {
    private:
        // M' typedef
        using MprimeMatrix = Eigen::Matrix<double, 6 - 3 * (4 - MESH_SIZE), 1>;

        // time-values
        std::array<double, MESH_SIZE> &values;

        // solver max-iteration and tolerance
        int max_iters;
        double tol;

        // M' matrix of the local solver
        MprimeMatrix M;

        // This method defines the M'-distance for Triangular elements
        inline double distance(const double l1) {
            if constexpr (MESH_SIZE == 3)
                return std::sqrt(l1 * l1 * M(0) + 2 * l1 * M(1) + M(2));
            else
                return 0;
        }

        // This method defines the M'-distance for Tetrahedral elements
        inline double distance(const double l1, const double l2) {
            if constexpr (MESH_SIZE == 4)
                return std::sqrt(l1 * l1 * M(0) + l2 * l2 * M(3) + 2 * (l1 * l2 * M(1) + l1 * M(2) + l2 * M(4)) + M(5));
            else
                return 0;
        }

        // This method perform the computation of M' matrix given velocity matrix V and element points
        inline void compute_MPrime(const Traits::VelocityM &V, const std::array<Traits::Point, MESH_SIZE> &points) {
            Eigen::Matrix<double, 3, MESH_SIZE - 1> E;

            E.col(0) << points[MESH_SIZE - 2] - points[0];
            E.col(1) << points[2] - points[1];

            if constexpr (MESH_SIZE == 4) {
                E.col(2) << points[3] - points[2];
                auto MI = E.transpose() * V * E;
                M << MI(0, 0), MI(0, 1), MI(0, 2), MI(1, 1), MI(1, 2), MI(2, 2);
            } else {
                auto MI = E.transpose() * (V * E);
                M << MI(0, 0), MI(0, 1), MI(1, 1);
            }
        }

    public:
        LocalSolver(const Traits::VelocityM &V, const std::array<Traits::Point, MESH_SIZE> &points,
                    std::array<double, MESH_SIZE> &values, int max_iters, double tol)
                : values(values), max_iters(max_iters), tol(tol) {
            compute_MPrime(V, points);
        }

        // This function performs the local solver, returning the solution value for the last point with reference to
        // point order of constructor input
        double operator()();
    };

    template
    class LocalSolver<3>;

    template
    class LocalSolver<4>;

}

#endif //LOCAL_SOLVER
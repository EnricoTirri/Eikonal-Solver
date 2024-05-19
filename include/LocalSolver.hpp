
#ifndef EIKONAL_SOLVER_LOCALSOLVER_H
#define EIKONAL_SOLVER_LOCALSOLVER_H

#include <Eigen/Core>
#include "EikonalTraits.hpp"

namespace Eikonal {

    template<int MESH_SIZE>
    class LocalSolver {
    private:
        using MprimeMatrix = Eigen::Matrix<double, 6 - 3 * (4 - MESH_SIZE), 1>;

        Eigen::Matrix<double, MESH_SIZE, 1> &values;
        int max_iters;
        double tol;

        MprimeMatrix M;

        inline double distance(const double l1) {
            if constexpr (MESH_SIZE == 3)
                return std::sqrt(l1 * l1 * M(0) + 2 * l1 * M(1) + M(2));
            else
                return 1;
        }

        inline double distance(const double l1, const double l2) {
            if constexpr (MESH_SIZE == 4)
                return std::sqrt(l1 * l1 * M(0) + l2 * l2 * M(3) + 2 * (l1 * l2 * M(1) + l1 * M(2) + l2 * M(4)) + M(5));
            else
                return 1;
        }

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
                    Eigen::Matrix<double, MESH_SIZE, 1> &values, int max_iters, double tol)
                : values(values), max_iters(max_iters), tol(tol) {
            compute_MPrime(V, points);
        }

        double operator()();
    };

    template
    class LocalSolver<3>;

    template
    class LocalSolver<4>;

}
#endif //EIKONAL_SOLVER_LOCALSOLVER_H

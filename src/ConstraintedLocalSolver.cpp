
#include "LocalSolver.hpp"
#include "EikonalMath.hpp"
#include <Eigen/LU>

namespace Eikonal {

    template<int MESH_SIZE>
    inline double LocalSolver<MESH_SIZE>::operator()() {
        if constexpr (MESH_SIZE == 3) {
            const double t12 = values[1] - values[0];
            SideFunction fun(M(0), M(1), M(2), t12);
            double l1 = constraintBisection(fun, max_iters, tol);
            values[2] = -t12 * l1 + values[1] + distance(l1);
            return values[2];
        } else if constexpr (MESH_SIZE == 4) {
            using Vector = Eigen::Matrix<double, 2, 1>;
            using Jacobian = Eigen::Matrix<double, 2, 2>;
            double l1 = 0;
            double l2 = 0;
            const double t13 = values[2] - values[0];
            const double t23 = values[2] - values[1];
            Vector R;
            int iters = max_iters;
            do {
                double dist = distance(l1, l2);
                double la = l1 * M(0) + l2 * M(1) + M(2);
                double lb = l1 * M(1) + l2 * M(3) + M(4);
                R << -t13 * dist + la, -t23 * dist + lb;
                Jacobian J;
                J.row(0) << M(0) - t13 * la / dist,
                        M(1) - t13 * lb / dist;
                J.row(1) << M(1) - t23 * la / dist,
                        M(3) - t23 * lb / dist;
                J = J.inverse().eval();
                Vector dir = -J * R;
                l1 += dir(0);
                l2 += dir(1);
                if (l1 <= 0) {
                    l1 = 0;
                    if (l2 <= 0) {
                        l2 = 0;
                    } else {
                        SideFunction fun(M(3), M(4), M(5), t23);
                        l2 = constraintBisection(fun, max_iters, tol);
                    }
                    break;
                } else {
                    if (l2 <= 0) {
                        SideFunction fun(M(0), M(2), M(5), t13);
                        l1 = constraintBisection(fun, max_iters, tol);
                        l2 = 0;
                        break;
                    } else if (l1 + l2 >= 1) {
                        double p1 = M(0) - 2 * M(1) + M(3);
                        double p2 = M(1) + M(2) - M(3) - M(4);
                        double p3 = M(3) + 2 * M(4) + M(5);
                        SideFunction fun(p1, p2, p3, t13 - t23);
                        l1 = constraintBisection(fun, max_iters, tol);
                        l2 = 1 - l1;
                        break;
                    }
                }
                iters--;
            } while (R.norm() > tol && iters != 0);
            values[3] = -t13 * l1 - t23 * l2 + values[2] + distance(l1, l2);
            return values[3];
        }
        return 0;
    }
}

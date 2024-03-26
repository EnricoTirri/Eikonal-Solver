//
// Created by Enrico on 23/03/2024.
//

#include "LocalSolver.hpp"
#include <Eigen/Core>
#include <Eigen/LU>

namespace Eikonal {

    class SideFunction {
    private:
        const double a, b, c, dt;

    public:
        SideFunction(const double a, const double b, const double c, const double dt) : a(a), b(b), c(c), dt(dt) {};

        inline double operator()(const double x) const {
            return (x * a + b) / std::sqrt(2 * x * b + x * x * a + c) - dt;
        }
    };

     int not_same_sign(const double a, double b){

        int sign_bit;
        asm (
                "xorpd %[a], %[b]\n\t"     // XOR the signs of a and b
                "vmovmskpd %[b], %[sign]\n\t"  // Move the sign bits to the lower 2 bits of b
                "and $1, %[sign]\n\t"      // Isolate the least significant bit
                : [b] "+x" (b), [sign] "=r" (sign_bit) //b is not really modified
        : [a] "x" (a)
        :
        );
        return sign_bit;
    }


    inline double constraintBisection(SideFunction f, const int max_iters, const double tol) {
        double a = 0;
        double b = 1;
        double fa = f(a);
        double fb = f(b);
        if (!not_same_sign(fa, fb)) {
            if (std::abs(fa) < std::abs(fb))
                return 0;
            else
                return 1;
        }
        double k;
        int iters = max_iters;
        double res;
        do {
            k = (a + b) / 2;
            res = f(k);
            if (!not_same_sign(fa, res)) {
                a = k;
                fa = res;
            } else {
                b = k;
            }
            iters--;
        } while (std::abs(res) > tol && iters != 0);
        return k;
    }

    template<int MESH_SIZE>
    double LocalSolver<MESH_SIZE>::operator()() {
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
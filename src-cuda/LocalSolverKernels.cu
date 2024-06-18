
#include "LocalSolverKernels.hpp"

namespace Eikonal {

    // EIKONAL MATH KERNELS --------------------------------------------------------------------------------
    __device__ inline bool sameSign(double a, double b) {
        return a * b >= 0;
    }

    __device__ inline double
    constraintBisection(double ma, double mb, double mc, double dt, const int max_iters, const double tol) {
        auto f = [ma, mb, mc, dt](const double x) {
            return (x * ma + mb) / std::sqrt(2 * x * mb + x * x * ma + mc) - dt;
        };

        double a = 0;
        double b = 1;
        double fa = f(a);
        double fb = f(b);
        if (sameSign(fa, fb)) {
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
            if (sameSign(fa, res)) {
                a = k;
                fa = res;
            } else {
                b = k;
            }
            iters--;
        } while (std::abs(res) > tol && iters != 0);
        return k;
    }
    // -------------------------------------------------------------------------------------------------------



    // GLOBALLY USED LOCAL SOLVER FUNCTIONS ----------------------------------------------------------
    __device__ int esgn(const int &kcode, const int &lcode, const int &scode) {
        return (2 * (scode < kcode) - 1) * (2 * (scode < lcode) - 1);
    }

    __device__ int signParity(const int &graycode) {
        int s = graycode - ((graycode >> 1) & 033333333333) - ((graycode >> 2) & 011111111111);
        return 1 - 2 * ((((s + (s >> 3)) & 030707070707) % 63) % 2);
    }
    // -----------------------------------------------------------------------------------------------



    // TRIANGULAR LOCAL SOLVER KERNELS --------------------------------------------------------------
    __device__ inline void getSigns3(int *signs, const int &shift) {
        signs[0] = signParity((24 >> shift) & 7);
        signs[1] = signParity((16 >> shift) & 7);
    }


    __device__ void getMprimeMatrix3(const int &ptidx, const MprimeMatrix<3> &MT, const double *valin,
                                     MprimeMatrix<3> &M, double *valout) {
        constexpr int RED_SIZE = 3 - 1;
        int gcodes[RED_SIZE];
        int signs[RED_SIZE];

        // get gcodes (rotated of -1)
        int refId = (3 + ptidx - 1) % 3;
        int refcode = (1 << refId);

        for (int i = 0; i < RED_SIZE; ++i) {
            int k = (ptidx + i) % 3;
            gcodes[(RED_SIZE + i - 1) % RED_SIZE] = refcode + (1 << k);
        }

        // get rotated vals (rotated of shift)
        int shift = RED_SIZE - ptidx;
        for (int i = 0; i < RED_SIZE; ++i) {
            valout[i] = valin[(3 + i - shift) % 3];
        }

        getSigns3(signs, shift);

        int e0 = gcodes[0] / 2 - 1;
        M(0) = MT(e0);

        int e1 = gcodes[1] / 2 - 1;
        M(1) = MT(e1);

        int g2 = gcodes[0] xor gcodes[1];
        int s01 = esgn(gcodes[0], gcodes[1], g2);
        int e2 = g2 / 2 - 1;
        M(2) = s01 * signs[0] * signs[1] * (MT(e0) + MT(e1) - MT(e2)) / 2;
    }

    __device__ inline double distance3(const double l1, const MprimeMatrix<3> &M) {
        return std::sqrt(l1 * l1 * M(0) + 2 * l1 * M(2) + M(1));
    }


    template<>
    __device__ double solveLocal<3>(const int &pointref, const MprimeMatrix<3> &MT, double *valin) {
        TTraits<3>::MprimeMatrix M;
        double valout[3];

        getMprimeMatrix3(pointref, MT, valin, M, valout);

        double t12 = valout[1] - valout[0];
        double l1 = constraintBisection(M(0), M(2), M(1), t12, 2, 10e-6);
        return -t12 * l1 + valout[1] + distance3(l1, M);
    }

    // ------------------------------------------------------------------------------------------------


    // TETRAHEDRAL LOCAL SOLVER KERNELS ---------------------------------------------------------------
    __device__ inline void getSigns4(int *signs, const int &shift) {
        signs[0] = signParity((80 >> shift) & 15);
        signs[1] = signParity((96 >> shift) & 15);
        signs[2] = signParity((192 >> shift) & 15);
    }

    __device__ inline double distance4(const double l1, const double l2, const MprimeMatrix<4> &M) {
        return std::sqrt(l1 * l1 * M(0) + l2 * l2 * M(1) + 2 * (l1 * l2 * M(2) + l1 * M(5) + l2 * M(4)) + M(3));
    }


    __device__ void getMprimeMatrix4(const int &ptidx, const MprimeMatrix<4> &MT, const double *valin,
                                     MprimeMatrix<4> &M, double *valout) {

        constexpr int RED_SIZE = 4 - 1;
        int gcodes[RED_SIZE];
        int signs[RED_SIZE];

        // get gcodes (rotated of -1)
        int refId = (4 + ptidx - 1) % 4;
        int refcode = (1 << refId);

        for (int i = 0; i < RED_SIZE; ++i) {
            int k = (ptidx + i) % 4;
            gcodes[(RED_SIZE + i - 1) % RED_SIZE] = refcode + (1 << k);
        }

        // get rotated vals (rotated of shift)
        int shift = RED_SIZE - ptidx;
        for (int i = 0; i < RED_SIZE; ++i) {
            valout[i] = valin[(4 + i - shift) % 4];
        }

        getSigns4(signs, shift);

        int e0 = gcodes[0] / 2 - 1;
        M(0) = MT(e0);

        int e1 = gcodes[1] / 2 - 1;
        M(1) = MT(e1);

        int g2 = gcodes[0] xor gcodes[1];
        int s01 = esgn(gcodes[0], gcodes[1], g2);
        int e2 = g2 / 2 - 1;
        M(2) = s01 * signs[0] * signs[1] * (MT(e0) + MT(e1) - MT(e2)) / 2;

        int e3 = gcodes[2] / 2 - 1;
        M(3) = MT(e3);

        int g4 = gcodes[1] xor gcodes[2];
        int s12 = esgn(gcodes[1], gcodes[2], g4);
        int e4 = g4 / 2 - 1;
        M(4) = s12 * signs[1] * signs[2] * (MT(e1) + MT(e3) - MT(e4)) / 2;

        int g5 = gcodes[0] xor gcodes[2];
        int s02 = esgn(gcodes[0], gcodes[2], g5);
        int e5 = g5 / 2 - 1;
        M(5) = s02 * signs[0] * signs[2] * (MT(e0) + MT(e3) - MT(e5)) / 2;
    }


    template<>
    __device__ double solveLocal<4>(const int &pointref, const MprimeMatrix<4> &MT, double *valin) {

        TTraits<4>::MprimeMatrix M;
        double valout[4];
        getMprimeMatrix4(pointref, MT, valin, M, valout);

        using Vector = Eigen::Matrix<double, 2, 1>;
        using Jacobian = Eigen::Matrix<double, 2, 2>;
        double l1 = 0;
        double l2 = 0;
        const double t13 = valout[2] - valout[0];
        const double t23 = valout[2] - valout[1];
        Vector R;
        int iters = 2;
        do {
            double dist = distance4(l1, l2, M);
            double la = l1 * M(0) + l2 * M(2) + M(5);
            double lb = l1 * M(2) + l2 * M(1) + M(4);
            R << -t13 * dist + la, -t23 * dist + lb;
            Jacobian J;
            J.row(0) << M(1) - t23 * lb / dist,
                        -(M(2) - t13 * lb / dist);
            J.row(1) << -(M(2) - t23 * la / dist),
                    M(0) - t13 * la / dist,

            J /= (J(0,0)*J(1,1) - J(0,1) * J(1,0));

            Vector dir = -J * R;
            l1 += dir(0);
            l2 += dir(1);
            if (l1 <= 0) {
                l1 = 0;
                if (l2 <= 0) {
                    l2 = 0;
                } else {
                    l2 = constraintBisection(M(1), M(4), M(3), t23, 5000, 10e-6);
                }
                break;
            } else {
                if (l2 <= 0) {
                    l1 = constraintBisection(M(0), M(5), M(3), t13, 5000, 10e-6);
                    l2 = 0;
                    break;
                } else if (l1 + l2 >= 1) {
                    double p1 = M(0) - 2 * M(2) + M(1);
                    double p2 = M(2) + M(5) - M(1) - M(4);
                    double p3 = M(1) + 2 * M(4) + M(3);
                    l1 = constraintBisection(p1, p2, p3, t13 - t23, 5000, 10e-6);
                    l2 = 1 - l1;
                    break;
                }
            }
            iters--;
        } while (R.norm() > 10e-6 && iters != 0);
        return -t13 * l1 - t23 * l2 + valout[2] + distance4(l1, l2, M);
    }



    // ------------------------------------------------------------------------------------------------
}

#include "LocalSolverKernels.hpp"

namespace Eikonal{

    // EIKONAL MATH KERNELS --------------------------------------------------------------------------------
    __device__ inline bool sameSign(double a, double b){
        return a*b >= 0;
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
    __device__ int esgn(const int &kcode, const int &lcode, const int &scode){
        return (2 * (scode < kcode) - 1) * (2 * (scode < lcode) - 1);
    }

    __device__ int signParity(const int &graycode){
        int s = s - ((s >> 1) & 033333333333) - ((s >> 2) & 011111111111);
        return 1 - 2 * ((((s + (s >> 3)) & 030707070707) % 63) % 2);
    }
    // -----------------------------------------------------------------------------------------------



    // TRIANGULAR LOCAL SOLVER KERNELS --------------------------------------------------------------
    __device__ inline void getSigns3(int *signs, const int &shift){
        signs[0] = signParity((24 >> shift) & 7);
        signs[1] = signParity((16 >> shift) & 7);
    }


    __device__ void getMprimeMatrix3(const int &ptidx, const MprimeMatrix<3> &MT, const double *valin,
                                        MprimeMatrix<3> &M, double *valout){
        constexpr int RED_SIZE = 3 - 1;
        int gcodes[RED_SIZE];
        int signs[RED_SIZE];

        // get gcodes (rotated of -1)
        int refId = (3 + ptidx - 1) % 3;
        int refcode = (1 << refId);

        for (int i = 0; i < RED_SIZE; ++i) {
            int k = (ptidx + i) % 3;
            gcodes[(RED_SIZE + i - 1)%RED_SIZE] = refcode + (1 << k);
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
    __device__ double solveLocal <3>(const int &pointref, const MprimeMatrix<3> &MT, double *valin){
        TTraits<3>::MprimeMatrix M;
        double valout[3];

        getMprimeMatrix3(pointref, MT, valin, M, valout);

        double t12 = valout[1] - valout[0];
        double l1 = constraintBisection(M(0), M(2), M(1), t12, 5000, 10e-6);
        return -t12 * l1 + valout[1] + distance3(l1, M);
    }

    // ------------------------------------------------------------------------------------------------


    // TETRAHEDRAL LOCAL SOLVER KERNELS ---------------------------------------------------------------
    template<>
    __device__ double solveLocal <4>(const int &pointref, const MprimeMatrix<4> &MT, double *valin){
//        TTraits<3>::MprimeMatrix M;
//        double valout[3];
//
//        getMprimeMatrix3(pointref, MT, valin, M, valout);
//
//        double t12 = valout[1] - valout[0];
//        double l1 = constraintBisection(M(0), M(2), M(1), t12, 5000, 10e-6);
//        return -t12 * l1 + valout[1] + distance3(l1, M);
        return 0;
    }



    // ------------------------------------------------------------------------------------------------
}
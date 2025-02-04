
#ifndef EIKONAL_MATH
#define EIKONAL_MATH

#include <cmath>
#ifdef __AVX2__
#include "immintrin.h"
#endif

namespace Eikonal {

    // This class represent the function of eikonal problem along the side of domain:
    // Lambda_1 + Lambda_2 <= 1,  0 <= Lambda_i <= 1
    class SideFunction {
    private:
        const double a, b, c, dt;

    public:
        SideFunction(const double a, const double b, const double c, const double dt) : a(a), b(b), c(c), dt(dt) {};

        inline double operator()(const double x) const {
            return (x * a + b) / std::sqrt(2 * x * b + x * x * a + c) - dt;
        }
    };

#ifdef __AVX2__
    inline int sameSign(double a, double b) {
        // Load the double 'a' and 'b' into 128-bit vectors
        __m128d va = _mm_set_sd(a);
        __m128d vb = _mm_set_sd(b);

        // Perform bitwise XOR on the vectors, this will set the sign bit if the signs are different
        __m128d vresult = _mm_xor_pd(va, vb);


        // Extract the sign bits of the result and not them, so that the sign bit is set if the signs are the same
        int sign = ~_mm_movemask_pd(vresult);

        // _mm_movemask_pd extract two sign bits, because it operates on packed doubles, we only need one bit
        return sign & 1;
    }
#else
    inline int sameSign(double a, double b) {
        return std::signbit(a) == std::signbit(b);
    }
#endif

    // This function performs a constrained bisection method on the function f 
    inline double constraintBisection(SideFunction f, const int max_iters, const double tol) {
        double a = 0;
        double b = 1;
        double fa = f(a);
        double fb = f(b);
        //if in the interval there is no intersection with zero axis, then
        //return the minimum value (is minimum value of derivative, so the lower the
        //value the closest the solution)
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

}

#endif //EIKONAL_MATH
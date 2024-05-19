
#ifndef EIKONAL_SOLVER_EIKONALMATH_H
#define EIKONAL_SOLVER_EIKONALMATH_H

#include <cmath>

namespace Eikonal{

    class SideFunction {
    private:
        const double a, b, c, dt;

    public:
        SideFunction(const double a, const double b, const double c, const double dt) : a(a), b(b), c(c), dt(dt) {};

        inline double operator()(const double x) const {
            return (x * a + b) / std::sqrt(2 * x * b + x * x * a + c) - dt;
        }
    };

    int not_same_sign(const double a, const double b) {
        int sign_bit;
        double temp;
        asm (
                "vpxor  %[a], %[b],%[temp]\n\t"     // XOR the signs of a and b
                "vmovmskpd %[temp], %[sign]\n\t"  // Move the sign bits to the lower 2 bits of b
                "and $1, %[sign]\n\t"      // Isolate the least significant bit
                :  [sign] "=r"(sign_bit),[temp] "=x"(temp)
        : [a] "x"(a), [b] "x"(b)
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

}

#endif //EIKONAL_SOLVER_EIKONALMATH_H

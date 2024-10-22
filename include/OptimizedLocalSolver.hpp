#ifndef OPTIMIZED_LOCAL_SOLVER
#define OPTIMIZED_LOCAL_SOLVER

#include "EikonalTraits.hpp"

namespace Eikonal {

    // This class defines the structure of the Eikonal Optimized Local solver
    template<int MESH_SIZE>
    class OptimizedLocalSolver {

        // M' typedef
        using MprimeMatrix = TTraits<MESH_SIZE>::MprimeMatrix;

        using MeshPoints = TTraits<MESH_SIZE>::MeshPoints;

        // solver max-iterations and tolerance
        int max_iters;
        double tol;

        // solver minimal-list of M' entries
        MprimeMatrix MT;

        // This method computes the M'-distance for Triangular elements
        inline double distance(const double l1, const MprimeMatrix &M) {
            if constexpr (MESH_SIZE == 3)
                return std::sqrt(l1 * l1 * M(0) + 2 * l1 * M(2) + M(1));
            else
                return 0;
        }

        // This methods computes the M'-distance for tetrahedral elements
        inline double distance(const double l1, const double l2, const MprimeMatrix &M) {
            if constexpr (MESH_SIZE == 4)
                return std::sqrt(l1 * l1 * M(0) + l2 * l2 * M(1) + 2 * (l1 * l2 * M(2) + l1 * M(5) + l2 * M(4)) + M(3));
            else
                return 0;
        }

        // This method computes the minimal list of M' entries
        void buildMprimeMatrix(const Traits::VelocityM &V, const MeshPoints &points) {
            Traits::Point E;
            E << points[1] - points[0];
            MT(0) = E.transpose() * V * E;

            E << points[2] - points[0];
            MT(1) = E.transpose() * V * E;

            E << points[2] - points[1];
            MT(2) = E.transpose() * V * E;

            if constexpr (MESH_SIZE == 4) {
                E << points[3] - points[0];
                MT(3) = E.transpose() * V * E;

                E << points[3] - points[1];
                MT(4) = E.transpose() * V * E;

                E << points[3] - points[2];
                MT(5) = E.transpose() * V * E;
            }
        }

        // Helper method to retrieve rotated M' configuration signs
        inline void getSigns(std::array<int, MESH_SIZE - 1> &signs, const int &shift) {
            if constexpr (MESH_SIZE == 3) {
                int s0 = (24 >> shift) & 7;
                s0 = s0 - ((s0 >> 1) & 033333333333) - ((s0 >> 2) & 011111111111);
                signs[0] = 1 - 2 * ((((s0 + (s0 >> 3)) & 030707070707) % 63) % 2);

                int s1 = (16 >> shift) & 7;
                s1 = s1 - ((s1 >> 1) & 033333333333) - ((s1 >> 2) & 011111111111);
                signs[1] = 1 - 2 * ((((s1 + (s1 >> 3)) & 030707070707) % 63) % 2);

            } else if constexpr (MESH_SIZE == 4) {
                int s0 = (80 >> shift) & 15;
                s0 = s0 - ((s0 >> 1) & 033333333333) - ((s0 >> 2) & 011111111111);
                signs[0] = 1 - 2 * ((((s0 + (s0 >> 3)) & 030707070707) % 63) % 2);

                int s1 = (96 >> shift) & 15;
                s1 = s1 - ((s1 >> 1) & 033333333333) - ((s1 >> 2) & 011111111111);
                signs[1] = 1 - 2 * ((((s1 + (s1 >> 3)) & 030707070707) % 63) % 2);
                // if constexpr (MESH_SIZE == 4) {
                int s2 = (192 >> shift) & 15;
                s2 = s2 - ((s2 >> 1) & 033333333333) - ((s2 >> 2) & 011111111111);
                signs[2] = 1 - 2 * ((((s2 + (s2 >> 3)) & 030707070707) % 63) % 2);
            }
        }

        // Helper method to retrieve mixed M' entries signs
        inline int esgn(const int &kcode, const int &lcode, const int &scode) {
            return (2 * (scode < kcode) - 1) * (2 * (scode < lcode) - 1);
        }


        // This method retrieve the rotated M' configuration with reference to a given ptIdx,
        // It also reorders the time-values with reference to the new configuration
        inline void
        getMprimeMatrix(const int &ptidx, const std::array<double, MESH_SIZE> &valin, MprimeMatrix &M, std::array<double, MESH_SIZE> &valout) {
            constexpr int RED_SIZE = MESH_SIZE - 1;
            std::array<int, RED_SIZE> gcodes;
            std::array<int, RED_SIZE> signs;

            // get gcodes (rotated of -1)
            int refId = (MESH_SIZE + ptidx - 1) % MESH_SIZE;
            int refcode = (1 << refId);

            for (int i = 0; i < RED_SIZE; ++i) {
                int k = (ptidx + i) % MESH_SIZE;
                gcodes[(RED_SIZE + i - 1)%RED_SIZE] = refcode + (1 << k);
            }

            // get rotated vals (rotated of shift)
            int shift = RED_SIZE - ptidx;
            for (int i = 0; i < RED_SIZE; ++i) {
                valout[i] = valin[(MESH_SIZE + i - shift) % MESH_SIZE];
            }

            getSigns(signs, shift);

            int e0 = gcodes[0] / 2 - 1;
            M(0) = MT(e0);

            int e1 = gcodes[1] / 2 - 1;
            M(1) = MT(e1);

            int g2 = gcodes[0] xor gcodes[1];
            int s01 = esgn(gcodes[0], gcodes[1], g2);
            int e2 = g2 / 2 - 1;
            M(2) = s01 * signs[0] * signs[1] * (MT(e0) + MT(e1) - MT(e2)) / 2;

            if(MESH_SIZE == 4){
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
        }


    public:
        OptimizedLocalSolver(int max_iters, double tol, const Traits::VelocityM &V, const MeshPoints &points)
                : max_iters(max_iters), tol(tol) {
            buildMprimeMatrix(V, points);
        }

        MprimeMatrix &getMprimeMatrix(){
            return MT;
        }

        // This function performs the local solver with reference to the given pointref,
        // returns the calculated value into the pointvalues list
        double operator()(const int &pointref, const std::array<double, MESH_SIZE> &pointValues);
    };

    template
    class OptimizedLocalSolver<3>;

    template
    class OptimizedLocalSolver<4>;
}

#endif //OPTIMIZED_LOCAL_SOLVER

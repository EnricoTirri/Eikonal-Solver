
#ifndef EIKONAL_TRAITS
#define EIKONAL_TRAITS

#include <Eigen/Core>

namespace Eikonal {

    // This class defines not-templated traits for eikonal global and local solvers
    class Traits {
    public:
        using Point = Eigen::Matrix<double, 3, 1>;
        using VelocityM = Eigen::Matrix<double, 3, 3>;
    };

    // This class defines templated traits for eikonal global and local solvers
    template<unsigned int MESH_SIZE>
    class TTraits{
    public:
        using MprimeMatrix = Eigen::Matrix<double, 6 - 3 * (4 - MESH_SIZE), 1>;
        using MeshPoints = std::array<Traits::Point, MESH_SIZE>;
    };
}

namespace std {

    // This struct defines the hash-operator for Eikonal Traits "Point"
    template<>
    struct hash<Eikonal::Traits::Point> {
        size_t operator()(const Eikonal::Traits::Point &k) const {
            size_t hashed;
            hashed = hash<double>()(k[0]);
            hashed ^= (hash<double>()(k[1]) + 0x9e3779b1 + (hashed << 6) + (hashed >> 2));
            hashed ^= (hash<double>()(k[2]) + 0x9e3779b1 + (hashed << 6) + (hashed >> 2));

            return hashed;
        }
    };
}

#endif //EIKONAL_TRAITS
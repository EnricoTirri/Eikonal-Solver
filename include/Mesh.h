//
// Created by giorgio on 20/11/23.
//

#ifndef EIKONEL_TEST_MESH_H
#define EIKONEL_TEST_MESH_H
#ifndef DIMENSION
#define DIMENSION 2
#define MESHSIZE 3
#else
#ifndef MESHSIZE
#define MESHSIZE 3
#endif
#endif

#include <Eigen/Core>
#include <vector>
#include <memory>

typedef Eigen::Matrix<double, DIMENSION, 1> Point;
typedef std::array<Point, MESHSIZE> MeshElement;

//now we will implement this algorithm Fast iterative method (X,L)
//define hash function for Point
namespace std {
    template<>
    struct hash<Point> {
        std::size_t operator()(const Point &k) const {
            using std::size_t;
            using std::hash;
            using std::string;
            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:
            std::size_t hashed = 0;
#if DIMENSION == 2
            hashed = hash<double>()(k[0]) ^ (hash<double>()(k[1]) + 0x9e3779b9 + (hashed << 6) + (hashed >> 2));
#else
            hashed = hash<double>()(k[0]) ^ (hash<double>()(k[1]) + 0x9e3779b9 + (hashed << 6) + (hashed >> 2))
       ^ (hash<double>()(k[2]) + 0x9e3779b9 + (hashed << 6) + (hashed >> 2));
#endif
            return hashed;
        }
    };

    template<>
    struct hash<MeshElement> {
        std::size_t operator()(const MeshElement &k) const {
            using std::size_t;
            using std::hash;
            using std::string;
            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:
            std::size_t hashed = 0;
#if MESHSIZE == 3
            hashed = hash<Point>()(k[0]) ^ (hash<Point>()(k[1]) + 0x9e3779b9 + (hashed << 6) + (hashed >> 2))
                     ^ (hash<Point>()(k[2]) + 0x9e3779b9 + (hashed << 6) + (hashed >> 2));
#else
            hashed = hash<Point>()(k[0]) ^ (hash<Point>()(k[1]) + 0x9e3779b9 + (hashed << 6) + (hashed >> 2))
       ^ (hash<Point>()(k[2]) + 0x9e3779b9 + (hashed << 6) + (hashed >> 2)) ^ (hash<Point>()(k[3]) + 0x9e3779b9 + (hashed << 6) + (hashed >> 2));
#endif
            return hashed;
        }
    };
}
typedef struct {
    std::size_t start, end;
} range_t;

class Mesh {
public:
    std::vector<MeshElement*> adjacentList;
    std::unordered_map<Point, range_t> index;
    std::vector<MeshElement> elements;
};


#endif //EIKONEL_TEST_MESH_H

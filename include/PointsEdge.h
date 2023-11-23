//
// Created by giorgio on 20/11/23.
//

#ifndef EIKONEL_TEST_POINTSEDGE_H
#define EIKONEL_TEST_POINTSEDGE_H


#include <Eigen/Core>
#include <vector>

typedef Eigen::Matrix<double, DIMENSION, 1> Point;
typedef std::array<Point, DIMENSION + 1> mesh_element;

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
            std::size_t hashed;
#if PHDIM == 2
            hashed = (hash<double>()(k[0])
                      ^ (hash<double>()(k[1]) << 1)) >> 1;
#else
            hashed = ((hash<double>()(k[0]) ^ (hash<double>()(k[1]) << 1)) ^ (hash<double>()(k[2]) << 2)>>2);
#endif
            return hashed;
        }
    };
}
typedef struct {
    std::size_t start, end;
} n_range;

class [[maybe_unused]] PointsEdge {
public:
    std::vector<mesh_element *> adjacentList;
    std::unordered_map<Point, n_range> index;
    std::vector<mesh_element> mesh;
};


#endif //EIKONEL_TEST_POINTSEDGE_H

//
// Created by giorgio on 20/11/23.
//

#ifndef EIKONEL_TEST_MESH_H
#define EIKONEL_TEST_MESH_H

#include <Eigen/Core>
#include <vector>
#include <memory>


using namespace Eikonal;


namespace std {
    template<>
    struct hash<Eikonal_traits<2u>::Point> {
        size_t operator()(const Eikonal_traits<2u>::Point &k) const {
            size_t hashed;
            hashed = hash<double>()(k[0]) ^ (hash<double>()(k[1]) + 0x9e3779b9 + (hashed << 6) + (hashed >> 2));

            return hashed;
        }
    };


    template<>
    struct hash<Eikonal_traits<3u>::Point> {
        size_t operator()(const Eikonal_traits<3u>::Point &k) const {
            size_t hashed;
            hashed = hash<double>()(k[0]);
            hashed ^= (hash<double>()(k[1]) + 0x9e3779b1 + (hashed << 6) + (hashed >> 2));
            hashed ^= (hash<double>()(k[2]) + 0x9e3779b1 + (hashed << 6) + (hashed >> 2));

            return hashed;
        }
    };
}

typedef struct {
    size_t start, end;
} range_t;

template<std::size_t DIM, std::size_t MESHSIZE>
class Mesh {
public:
    using MeshElement = std::array<int,MESHSIZE>;

    std::vector<int> adjacentList;
    std::vector<range_t> index;
    std::vector<typename Eikonal_traits<DIM>::Point> points;
    std::vector<MeshElement> elements;
};


#endif

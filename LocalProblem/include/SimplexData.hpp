/*
 * SimplexData.hpp
 *
 *  Created on: Jun 17, 2022
 *      Author: forma
 */

#ifndef EIKONAL_SIMPLEXDATA_HPP_
#define EIKONAL_SIMPLEXDATA_HPP_

#include "Eikonal_traits.hpp"

namespace Eikonal {
//#define PHDIM 3
//#define MSHDIM 3
    template<std::size_t PHDIM, std::size_t MSHDIM>
    struct SimplexData {
        using AnisotropyM = typename Eikonal_traits<PHDIM>::AnisotropyM;
        using MMatrix = typename Eikonal_traits<MSHDIM - 1>::MMatrix;
        using EMatrix = typename Eigen::Matrix<double,PHDIM,MSHDIM-1>;
        using Point = typename Eikonal_traits<PHDIM>::Point;


        std::array<Point, MSHDIM> points;
        MMatrix MM_Matrix;
        EMatrix E;

        //! This constructor just takes of vector that describes the simplex

        SimplexData(std::array<std::array<double, PHDIM>, MSHDIM> const &p,
                    AnisotropyM const &M = AnisotropyM::Identity()) {
            for (auto i = 0u; i < MSHDIM; ++i)
                points[i] = Eigen::Map<Point>(const_cast<double *>(p[i].data()));
            setup(M);
        }

        SimplexData(std::array<Point, MSHDIM> const &p,
                    AnisotropyM const &M = AnisotropyM::Identity()) :
                points{p} {
            setup(M);
        };

    private:
        void setup(AnisotropyM const &M) {
            E.col(0) = points[MSHDIM - 2u] - points[0];//e13 or e12
            if constexpr (MSHDIM == 4u) {
                E.col(1) = points[MSHDIM - 2u] - points[1]; //e23
            }
            E.col(MSHDIM - 2) = points[MSHDIM - 1u] - points[MSHDIM - 2u]; //e34 or e23
            MM_Matrix = E.transpose() * M * E;
        }
    };
}


#endif /* EIKONAL_SIMPLEXDATA_HPP_ */

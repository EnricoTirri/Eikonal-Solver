/*
 * Phi.hpp
 *
 *  Created on: Jun 16, 2022
 *      Author: forma
 */

#ifndef EIKONAL_PHI_HPP_
#define EIKONAL_PHI_HPP_

#include "LineSearch_traits.hpp"
#include "SimplexData.hpp"
#include <cmath>

namespace Eikonal {
    /**
     * @brief Phi function for Local Problem solution
     *
     * @tparam PHDIM the physical dimension
     * @tparam MSHDIM the mesh number of points
     */
//#define MSHDIM 3
//#define MSHDIM 3

    template<std::size_t PHDIM, std::size_t MSHDIM>
    struct Phi {
        using Scalar = double;
        // The dimension of the unknowns (Problem dimension - 2)
        static constexpr std::size_t DIM = MSHDIM - 2u;

    public:
        //@todo implementare move semantic per velocizzare
        Phi(SimplexData<PHDIM, MSHDIM> const &simplex, Eikonal_traits<MSHDIM>::Vector const &values) : simplexData{simplex}, values{values} {
            du(0) = values(0) - values(DIM);// u31 (u21)
            du(DIM) = values(DIM); // u3 (u2)
            if constexpr (MSHDIM == 4u) {
                du(1) = values(1) - values(DIM);
            }
            lambdaExt(DIM) = 1.0;
        };

        // a more efficient implementation chaches some quantities
        // here I just use the expressions straight away to avoid
        // errors
        Scalar operator()(Eikonal_traits<DIM>::Vector const &v) const {
            lambdaExt.template topRows<DIM>() = v;
            return lambdaExt.dot(du) + normL();
        }

        decltype(auto) gradient(Eikonal_traits<DIM>::Vector const &v) const {
            lambdaExt.template topRows<DIM>() = v;
            return du.template topRows<DIM>() +
                   simplexData.MM_Matrix.template block<DIM, DIM+1>(0, 0) * lambdaExt / normL();
        }

        typename Eikonal_traits<DIM>::MMatrix hessian(typename Eikonal_traits<DIM>::Vector const &v) const {
            lambdaExt.template topRows<DIM>() = v;
            auto n = 1. / normL();
            auto n3 = n * n * n;
            typename Eikonal_traits<DIM>::Vector part{simplexData.MM_Matrix.template block<DIM, DIM+1>(0, 0) * lambdaExt};
            typename Eikonal_traits<DIM>::MMatrix parta{simplexData.MM_Matrix.template block<DIM, DIM>(0, 0)};
            typename Eikonal_traits<DIM>::MMatrix partb{part * part.transpose()};
            return n * parta - n3 * partb;

        }

        Scalar normL() const {
            return std::sqrt(
                    lambdaExt.transpose() *
                    simplexData.MM_Matrix * lambdaExt
            );
        }

    private:
        SimplexData<PHDIM, MSHDIM> simplexData;
        typename Eikonal_traits<MSHDIM>::Vector values;
        typename Eikonal_traits<MSHDIM - 1>::Vector du;

        mutable typename Eikonal_traits<
                MSHDIM - 1>::Vector lambdaExt;    // can change value from a const function because declared as mutable
    };
}


#endif /* EIKONAL_PHI_HPP_ */

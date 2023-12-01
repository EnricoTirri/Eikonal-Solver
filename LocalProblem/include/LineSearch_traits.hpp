/*
 * LineSearch_traits.hpp
 *
 *  Created on: Dec 27, 2020
 *      Author: forma
 */

#ifndef EXAMPLES_SRC_LINESEARCH_LINESEARCH_TRAITS_HPP_
#define EXAMPLES_SRC_LINESEARCH_LINESEARCH_TRAITS_HPP_

#include "Eigen/Core"
#include <functional>

namespace apsc::LineSearch {

    using Scalar = double;

    template<std::size_t DIM>
    struct Traits {
        using Vector = Eigen::Matrix<Scalar, DIM, 1>;
        using Matrix = Eigen::Matrix<Scalar, DIM, DIM>;
        using CostFunction = std::function<Scalar(Vector const &)>;
        using Gradient = std::function<Vector(Vector const &)>;
        using Hessian = std::function<Matrix(Vector const &)>;
    };
}


#endif /* EXAMPLES_SRC_LINESEARCH_LINESEARCH_TRAITS_HPP_ */

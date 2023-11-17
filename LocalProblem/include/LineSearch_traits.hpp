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
namespace apsc
{
template<std::size_t DIM>
struct LineSearch_traits_base
{
  using Scalar = double;
  using Vector = Eigen::Matrix<Scalar, DIM, 1>;
  using Matrix = Eigen::Matrix<Scalar, DIM, DIM>;
  using CostFunction = std::function<Scalar(Vector const &)>;
  using Gradient = std::function<Vector(Vector const &)>;
  using Hessian  = std::function<Matrix(Vector const &)>;
};

#if DIMENSION==2
//#warning "two dimensional problem"
using LineSearch_traits=LineSearch_traits_base<1u>;     //solver traits base con DIM = 1, 1u sta per unsigned int
#else
using LineSearch_traits=LineSearch_traits_base<2u>;
#endif

/*
struct LineSearch_traits
{
  using Scalar = double;
  using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using CostFunction = std::function<Scalar(Vector const &)>;
  using Gradient = std::function<Vector(Vector const &)>;
  using Hessian  = std::function<Matrix(Vector const &)>;
};
*/
} // namespace apsc

#endif /* EXAMPLES_SRC_LINESEARCH_LINESEARCH_TRAITS_HPP_ */

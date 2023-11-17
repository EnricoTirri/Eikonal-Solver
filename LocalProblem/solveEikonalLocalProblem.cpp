/*
 * solveEikonalLocalProblem.cpp
 *
 *  Created on: Jun 18, 2022
 *      Author: forma
 */
#include "solveEikonalLocalProblem.hpp"
#if DIMENSION==2
template class Eikonal::solveEikonalLocalProblem<2u>;
constexpr std::size_t PHDIM=2u;
#else
template class Eikonal::solveEikonalLocalProblem<3u>;
constexpr std::size_t PHDIM=3u;
#endif
//! The namespace of pyEikonal sofware
namespace Eikonal
{
  /*!
    Solves the local optimization problem using a line search algorithm
    @tparam PHDIM The space dimension (2 or 3) Given here as cpp Macro (@todo can be changed to a tempalte parameter)
    @param  the coordinated of the vertices of a  a triangular/tetrahedral element ordered so that the last vertex 
            is the one where we have to compute u
    @param values. An Eigen vector containing the values of u at the base of the tiangle/tetra
    @param M A positive definite matrix (as an Eigen Matrix) the contains the information of the wave celerity
   */
 EikonalSolution<PHDIM> solveLocalProblem(std::array<std::array<double,PHDIM>,PHDIM+1u> element,
		 Eigen::Matrix<double,PHDIM,1>values,
         Eigen::Matrix<double,PHDIM,PHDIM> const & M)
 {
	 Eikonal::SimplexData<PHDIM> simplex{element,M};
	 Eikonal::solveEikonalLocalProblem<PHDIM> solver{std::move(simplex),std::move(values)};
	 return solver();
 }
  /*!
`@brief The version for the isotropica case with celerity=1
   */
 EikonalSolution<PHDIM> solveLocalProblemIsotropic(std::array<std::array<double,PHDIM>,PHDIM+1u> element,
 		 Eigen::Matrix<double,PHDIM,1> values)
 {
	 return solveLocalProblem(element,values,Eigen::Matrix<double,PHDIM,PHDIM>::Identity());
 }
}
//template class Eikonal::solveEikonalLocalProblem<2u>;




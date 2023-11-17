#ifndef HH_SOLVEEIKONAL__HH
#define HH_SOLVEEIKONAL__HH
#include "DescentDirectionFactory.hpp"
#include "DescentDirections.hpp"
#include "GradientFiniteDifference.hpp"
#include "LineSearch.hpp"
#include "LineSearchSolver.hpp"
#include "Phi.hpp"
#include <cmath>
#include <iostream>
#include <memory>
namespace Eikonal
{
/*!
 * @brief Structure returning the solution of the local eikonal problem
 *
 * @tparam PHDIM The physical dimension
 */
template <std::size_t PHDIM> struct EikonalSolution
{
  using Vector=typename apsc::LineSearch_traits_base<PHDIM - 1u>::Vector;
  double  value; //! The value at the new point
  Vector  lambda; //! The value(s) of lambda (foot of the characteristics)
  int     status; //! 0= converged 1=no descent direction 2=no convergence
};

/*!
 * @brief Driver for the local solver
 *
 * @tparam PHDIM The physical dimension
 */
template <std::size_t PHDIM> class solveEikonalLocalProblem
{
public:
  using Vector = apsc::LineSearch_traits::Vector;
  using Matrix = apsc::LineSearch_traits::Matrix;
  //! I pass a simplex structure and the values in the constructor
  //! @todo To save memory and time I have to store references in Phi
  template <typename SIMPLEX, typename VALUES>
  solveEikonalLocalProblem(SIMPLEX &&simplex, VALUES &&values)
    : my_phi{std::forward<SIMPLEX>(simplex), std::forward<VALUES>(values)}
  {}
  /*!
   * Solves the local problem
   */
  EikonalSolution<PHDIM>          //restituisce una struct
  operator()() const
  {
    apsc::OptimizationData optimizationData;
    optimizationData.NumberOfVariables = PHDIM - 1u;
    optimizationData.costFunction = [this](const Vector &x) {
      return this->my_phi(x);
    };
    optimizationData.gradient = [this](const Vector &x) {
      return this->my_phi.gradient(x);
    };

    optimizationData.hessian = [this](const Vector &x) {
      return this->my_phi.hessian(x);
    };
    Vector initialPoint;
    if constexpr(PHDIM == 3u)
      {
        setBounds(optimizationData, {0., 0.}, {1.0, 1.0});
      }
    else
      {
        setBounds(optimizationData, {0.0}, {1.0});
      }
    initialPoint.fill(0.333);           //put lambdas = 1/3 all elements in the vector, it is the initial value
    apsc::LinearSearchSolver solver(optimizationData,
                                    std::make_unique<apsc::NewtonDirection>(),
                                    optimizationOptions, lineSearchOptions);
    solver.setInitialPoint(initialPoint);
    auto [finalValues, numIter, status] = solver.solve();
#ifdef VERBOSE
    if(status == 0)
      std::cout << "Solver converged" << std::endl;
    else
      std::cout << "Solver DID NOT converge" << std::endl;

    std::cout << "Point found=" << finalValues.currentPoint.transpose()
              << "\nCost function value=" << finalValues.currentCostValue
              << "\nGradient   =" << finalValues.currentGradient.transpose()
              << "\nGradient norm=" << finalValues.currentGradient.norm()
              << "\nNumber of iterations=" << numIter << "\nStatus=" << status
              << std::endl;
#endif
    return {finalValues.currentCostValue, finalValues.currentPoint, status};
  }
  static void setLineSearchOptions(apsc::LineSearchOptions const & lso)
		{
	  lineSearchOptions=lso;
		}
  static void setOptimizationOptions(apsc::OptimizationOptions const & oop)
 		{
 	  optimizationOptions=oop;
 		}
private:
  Eikonal::Phi<PHDIM>                     my_phi;
  inline static apsc::LineSearchOptions   lineSearchOptions;
  inline static apsc::OptimizationOptions optimizationOptions;
  inline static std::unique_ptr<apsc::DescentDirectionBase>
    descentDirectionFunction = std::make_unique<apsc::NewtonDirection>();
};
#if DIMENSION==2
extern template class solveEikonalLocalProblem<2u>;
/*!
 * @brief A function that calls the solution of the local problem
 *
 * @param element An array \f$d+1\times d\f$ that contains the points of the simplex (tetra or triangle) with the convention that
 * the first three rows contain the coordinate of the base, tha last row that of the forth vertex where the solution is unknown
 * @param values The values of t at the base points
 * @param M The anisotropy matrix.
 * @return A structure containing the solution
 */
EikonalSolution<2u> solveLocalProblem(std::array<std::array<double,2u>,3u> element,
		Eigen::Matrix<double,2u,1u> values,
         Eigen::Matrix<double,2u,2u> const & M);
/*!
 * @brief The version of solveLocalProblem for the standard probelm where M=I (identity)
 *
 * @detail See the documantation of the other version
 */
EikonalSolution<2u> solveLocalProblemIsotropic(std::array<std::array<double,2u>,3u> element,
		 Eigen::Matrix<double,2u,1u> values);
#else
extern template class solveEikonalLocalProblem<3u>;
/*!
 * @brief A function that calls the solution of the local problem
 *
 * @param element An array \f$d+1\times d\f$ that contains the points of the simplex (tetra or triangle) with the convention that
 * the first three rows contain the coordinate of the base, tha last row that of the forth vertex where the solution is unknown
 * @param values The values of t at the base points
 * @param M The anisotropy matrix.
 * @return A structure containing the solution
 */
EikonalSolution<3u> solveLocalProblem(std::array<std::array<double,3u>,4u> element,
		 Eigen::Matrix<double,3u,1u> values,
         Eigen::Matrix<double,3u,3u> const & M);
/*!
 * @brief The version of solveLocalProblem for the standard probelm where M=I (identity)
 *
 * @detail See the documantation of the other version
 */
EikonalSolution<3u> solveLocalProblemIsotropic(std::array<std::array<double,3u>,4u> element,
		Eigen::Matrix<double,3u,1u> values);
#endif

} // namespace Eikonal
#endif

/*
 * DescentDirections.hpp
 *
 *  Created on: Dec 28, 2020
 *      Author: forma
 */

#ifndef EXAMPLES_SRC_LINESEARCH_DESCENTDIRECTIONS_HPP_
#define EXAMPLES_SRC_LINESEARCH_DESCENTDIRECTIONS_HPP_
#include "DescentDirectionBase.hpp"
#include "Eigen/Dense"
#include <limits>
#include <numeric>
#include <functional>
#include <array>
namespace apsc
{
/*!
 * Implements the gradient search
 */
class GradientDirection : public DescentDirectionBase
{
public:
  /*!
   *  Returns - gradient as descent direction
   * @param values The current values.
   * @return The descent direction.
   */
  apsc::LineSearch_traits::Vector
  operator()(apsc::OptimizationCurrentValues const &values) override
  {
    return -values.currentGradient;
  }
  /*!
   * @brief The class is clonable
   *
   * @return A clone of myself wrapped into a unique pointer
   */
  virtual
  std::unique_ptr<DescentDirectionBase>
  clone() const override
  {return std::make_unique<GradientDirection>(*this);}

};
/*!
 * Implements the gradient search
 */
class NewtonDirection : public DescentDirectionBase
{
public:
	/*!
	 *  Returns - gradient as descent direction
	 * @param values The current values.
	 * @return The descent direction.
	 */
	apsc::LineSearch_traits::Vector
	operator()(apsc::OptimizationCurrentValues const &values) override
	{
		//    if(!values.bounded)
		//      return values.currentHessian.llt().solve(-values.currentGradient);
		//    else
		//      {
#if DIMENSION == 3
		bool active=false;
		std::array<bool,3> constrained={false,false,false};
		for (auto i=0; i<values.currentPoint.size();++i)
		{
			constrained[i]=
					(values.currentPoint[i]==0.0
							and values.currentGradient[i]>0 );
			active= active || constrained[i];
		}

		constrained[2]=(std::abs(values.currentPoint[0]+values.currentPoint[1]-1)<=eps)
				and (values.currentGradient[0]+values.currentGradient[1])<=0.0;
		active=active||constrained[2];
#else
		bool active=((values.currentPoint[0]==0.0
				and values.currentGradient[0]>0) or
				(values.currentPoint[0]==1.0
				and values.currentGradient[0]<0));
		if(active)
			return ZeroVec;
		else
			return -values.currentHessian.inverse()*values.currentGradient;
#endif
#if DIMENSION == 3
		if(not active)
			return -values.currentHessian.inverse()*values.currentGradient;
		if(
				(constrained[0] and constrained[1])
				or
				(values.currentPoint[0]==1. and values.currentPoint[1]==0. and
						(values.currentGradient[1]-values.currentGradient[0])>=0. and
						values.currentGradient[0]<=0.)
						or
						(values.currentPoint[0]==0. and values.currentPoint[1]==1. and
								(values.currentGradient[1]-values.currentGradient[0])<=0. and
								values.currentGradient[1]<=0.)
		)
		{
			// gradient is pushing outside the constrained area
			return ZeroVec;
		}

		apsc::LineSearch_traits::Matrix Hi =values.currentHessian.inverse();
		if(constrained[0])
		{
			Hi.row(0).fill(0.);
			Hi.col(0).fill(0.);
		}
		else if	(constrained[1])
		{
			Hi.row(1).fill(0.);
			Hi.col(1).fill(0.);
		}
		else if(constrained[2])
		{
			Hi=P3*Hi*P3;
		}
	return -Hi*values.currentGradient;
#endif
}

/*!
 * @brief The class is clonable
 *
 * @return A clone of myself wrapped into a unique pointer
 */
virtual
std::unique_ptr<DescentDirectionBase>
clone() const override {return std::make_unique<NewtonDirection>(*this);}

private:
static inline const apsc::LineSearch_traits::Vector ZeroVec{apsc::LineSearch_traits::Vector::Zero()};
static inline const apsc::LineSearch_traits::Matrix P3{apsc::LineSearch_traits::Matrix::Identity()-
	0.5*apsc::LineSearch_traits::Matrix::Ones()};
static constexpr double eps=100.*std::numeric_limits<double>::epsilon();
};
/*!
 *  Implements the classic BFGS quasi-Newton algorithm.
 */
class BFGSDirection : public DescentDirectionBase
{
public:
  apsc::LineSearch_traits::Vector
  operator()(apsc::OptimizationCurrentValues const &values) override;
  /*!
   * You need to reset if you run another problem or you start from a different
   * initial point.
   * @note This is done inside LineSearchSolver
   */
  void reset() override;
  /*!
    * @brief The class is clonable
    *
    * @return A clone of myself wrapped into a unique pointer
    */
   virtual
   std::unique_ptr<DescentDirectionBase>
   clone() const override
   {return std::make_unique<BFGSDirection>(*this);}


private:
  apsc::OptimizationCurrentValues previousValues;
  Eigen::MatrixXd                 H;
  bool                            firstTime{true};
  double const smallNumber = std::sqrt(std::numeric_limits<double>::epsilon());
};

/*!
 * Implements BFGS with the direct computation of the approximate inverse of
 * the Hessian.
 */
class BFGSIDirection : public DescentDirectionBase
{
public:
  apsc::LineSearch_traits::Vector
       operator()(apsc::OptimizationCurrentValues const &values) override;
  void reset() override;

private:
  apsc::OptimizationCurrentValues previousValues;
  Eigen::MatrixXd                 H;
  bool                            firstTime{true};
  double const smallNumber = std::sqrt(std::numeric_limits<double>::epsilon());

  /*!
    * @brief The class is clonable
    *
    * @return A clone of myself wrapped into a unique pointer
    */
   virtual
   std::unique_ptr<DescentDirectionBase>
   clone() const override
   {return std::make_unique<BFGSIDirection>(*this);}


};
/*!
 * Bazrzilain-Borwein
 */
class BBDirection : public DescentDirectionBase
{
public:
  apsc::LineSearch_traits::Vector
  operator()(apsc::OptimizationCurrentValues const &values) override;
  void
  reset() override
  {
    firstTime = true;
  };
  /*!
    * @brief The class is clonable
    *
    * @return A clone of myself wrapped into a unique pointer
    */
   virtual
   std::unique_ptr<DescentDirectionBase>
   clone() const override
   {return std::make_unique<BBDirection>(*this);}


private:
  apsc::OptimizationCurrentValues previousValues;
  bool                            firstTime{true};
  double const smallNumber = std::sqrt(std::numeric_limits<double>::epsilon());
};
/*!
 * Non-linear CG with Polak-Ribiere formula
 */
class CGDirection : public DescentDirectionBase
{
public:
  apsc::LineSearch_traits::Vector
  operator()(apsc::OptimizationCurrentValues const &values) override;
  void
  reset() override
  {
    firstTime = true;
  };

  /*!
    * @brief The class is clonable
    *
    * @return A clone of myself wrapped into a unique pointer
    */
   virtual
   std::unique_ptr<DescentDirectionBase>
   clone() const override
   {return std::make_unique<CGDirection>(*this);}


private:
  apsc::OptimizationCurrentValues previousValues;
  bool                            firstTime{true};
  //! I need to keep track of previous descent direction
  apsc::LineSearch_traits::Vector prevDk;
};

} // namespace apsc

#endif /* EXAMPLES_SRC_LINESEARCH_DESCENTDIRECTIONS_HPP_ */

/*
 * Optimization_options.hpp
 *
 *  Created on: Dec 27, 2020
 *      Author: forma
 */
#ifndef EXAMPLES_SRC_LINESEARCH_OPTIMIZATION_OPTIONS_HPP_
#define EXAMPLES_SRC_LINESEARCH_OPTIMIZATION_OPTIONS_HPP_

#include "LineSearch_traits.hpp"
#include <exception>
#include <vector>

namespace apsc {
/*!
 * It holds all options for the line search algorithm
 */
    struct OptimizationOptions {
        using Scalar = apsc::LineSearch::Scalar;
        Scalar relTol = 1.e-5; //!< relative tolerance
        Scalar absTol = 1.e-5; //!< absolute tolerance
        unsigned int maxIter = 500;  //!< max n. of Iteration
        //@todo method to read from file
    };

/*!
 * It holds the main data for the optimization algorithm
 */
    template<size_t PRDIM>
    struct OptimizationData {
        apsc::LineSearch::Traits<PRDIM>::CostFunction costFunction; //!< The cost function.
        apsc::LineSearch::Traits<PRDIM>::Gradient gradient;              //!< The gradient of the cost function.
        //! The Hessian: by default an empty matrix
        apsc::LineSearch::Traits<PRDIM>::Hessian hessian = [](apsc::LineSearch::Traits<PRDIM>::Vector const &) {
            return typename apsc::LineSearch::Traits<PRDIM>::Matrix{};
        };
        //! The number of variables of the problem.
        std::size_t NumberOfVariables = 0;
        //! We may have bound constraints
        bool bounded = false;
        //! Uses if bounded=true
        std::vector<double> lowerBounds;
        //! Uses if bounded=true
        std::vector<double> upperBounds;

        static void setBounds(OptimizationData<PRDIM> &optimizationData, std::vector<double> const &lo,
                       std::vector<double> const &up) {
            if (optimizationData.NumberOfVariables > lo.size() and
                optimizationData.NumberOfVariables > up.size()) {
                throw std::runtime_error("Wrong bound sizes");
            }
            optimizationData.bounded = true;
            optimizationData.lowerBounds = lo;
            optimizationData.upperBounds = up;
        }
    };

/*!
 * A structure used to hold the current values.
 */
    template<size_t PRDIM>
    struct OptimizationCurrentValues {
        apsc::LineSearch::Scalar currentCostValue; //!< current cost.
        apsc::LineSearch::Traits<PRDIM>::Vector currentPoint;     //!< current point.
        apsc::LineSearch::Traits<PRDIM>::Vector currentGradient;  //!< current gradient.
        apsc::LineSearch::Traits<PRDIM>::Matrix currentHessian;  //!< current Hessian.
        bool bounded = false;       //! We may have bound constraints
        std::vector<double> lowerBounds;   //! Uses if bounded=true
        std::vector<double> upperBounds;   //! Uses if bounded=true
    };
} // namespace apsc

#endif /* EXAMPLES_SRC_LINESEARCH_OPTIMIZATION_OPTIONS_HPP_ */

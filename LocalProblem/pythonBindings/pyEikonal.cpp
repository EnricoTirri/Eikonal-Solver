/*
 * pyEikonal.cpp
 *
 *  Created on: Jun 22, 2022
 *      Author: forma
 */
#include "LineSearch_options.hpp"
#include "Optimization_options.hpp"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/eigen.h"
#include "solveEikonalLocalProblem.hpp"
using namespace apsc;
namespace py = pybind11;
#if DIMENSION == 3
static constexpr std::size_t PHDIM = 3u;
#define MODULENAME pyEikonal3
#else
static constexpr std::size_t PHDIM = 2u;
#define MODULENAME pyEikonal2
#endif
using Matrix = Eigen::Matrix<double, PHDIM, PHDIM>;
using Element = std::array<std::array<double, PHDIM>, PHDIM + 1u>;
using Values = std::array<double, PHDIM>;
using Solution = Eikonal::EikonalSolution<PHDIM>;
using Lambda = apsc::LineSearch_traits_base<PHDIM - 1u>::Vector;
using Solver = Eikonal::solveEikonalLocalProblem<PHDIM>;
/*!
 * @brief Binding module to python
 */
PYBIND11_MODULE(MODULENAME, m)
{
  py::class_<LineSearchOptions>(m, "LineSearchOptions")
    .def(py::init<>())
    .def_readwrite("sufficientDecreaseCoefficient",
                   &LineSearchOptions::sufficientDecreaseCoefficient,
				   "Coefficient for the sufficient decrease condition")
    .def_readwrite("stepSizeDecrementFactor",
                   &LineSearchOptions::stepSizeDecrementFactor,
				   "Factor of reduction of the step during backtracing")
    .def_readwrite("secondWolfConditionFactor",
                   &LineSearchOptions::secondWolfConditionFactor,"Not used")
    .def_readwrite("initialStep", &LineSearchOptions::initialStep,
    		"The initial step length (should be 1.0)")
    .def_readwrite("maxIter", &LineSearchOptions::maxIter,
    		"The maximum number of backtracing steps");

  py::class_<OptimizationOptions>(m, "OptimizationOptions")
    .def(py::init<>())
    .def_readwrite("relTol", &OptimizationOptions::relTol,
    		"Relative tolerance, used for the condition on the gradient")
    .def_readwrite("absTol", &OptimizationOptions::absTol,
    		"Absolute tolerance. Used for the condition on the step")
    .def_readwrite("maxIter", &OptimizationOptions::maxIter,
    		"Max number of iteration of the line search algorithm");

  py::class_<Solution>(m, "EikonalSolution")
    .def(py::init<>())
    .def_readwrite("value", &Solution::value,"The found value")
    .def_readwrite("l", &Solution::lambda,
    		"The lambda coordinate of the foot of the characteristic")
    .def_readwrite("status", &Solution::status,"Status: 0 means ok");

  m.def("setLineSearchOptions", [](LineSearchOptions const &opt) {
    Solver::setLineSearchOptions(opt);
  },"Sets line search options");

  m.def("setOptimizationOptions", [](OptimizationOptions const &opt) {
    Solver::setOptimizationOptions(opt);
  },"Sets general optimization options");

  m.def("solveLocalProblem", &Eikonal::solveLocalProblem,
        "A function that solves the local Eikonal problem", py::arg("element"),
        py::arg("values"), py::arg("M"));

  m.def("solveLocalProblemIsotropic", &Eikonal::solveLocalProblemIsotropic,
        "A function that solves the local Eikonal problem (isotropic)",
        py::arg("element"), py::arg("values"));
}

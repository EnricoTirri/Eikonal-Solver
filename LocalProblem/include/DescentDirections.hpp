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

namespace apsc {

    template<size_t PRDIM>
    struct DescentDirection {
        using Vector = apsc::LineSearch::Traits<PRDIM>::Vector;
        using OptimizationCurrentValues = apsc::OptimizationCurrentValues<PRDIM>;
        using Matrix = apsc::LineSearch::Traits<PRDIM>::Matrix;

/*!
 * Implements the gradient search
 */
        class GradientDirection : public DescentDirectionBase<PRDIM> {
        public:
            /*!
             *  Returns - gradient as descent direction
             * @param values The current values.
             * @return The descent direction.
             */
            Vector
            operator()(OptimizationCurrentValues const &values) override {
                return -values.currentGradient;
            }

            /*!
             * @brief The class is clonable
             *
             * @return A clone of myself wrapped into a unique pointer
             */
            virtual
            std::unique_ptr<DescentDirectionBase<PRDIM>>
            clone() const override { return std::make_unique<GradientDirection>(*this); }

        };

/*!
 * Implements the gradient search
 */

        class NewtonDirection : public DescentDirectionBase<PRDIM> {
        public:
            /*!
             *  Returns - gradient as descent direction
             * @param values The current values.
             * @return The descent direction.
             */
            Vector
            operator()(OptimizationCurrentValues const &values) override {
                //    if(!values.bounded)
                //      return values.currentHessian.llt().solve(-values.currentGradient);
                //    else
                //      {
                if constexpr (PRDIM == 2) {
                    bool active = false;
                    std::array<bool, 3> constrained = {false, false, false};
                    for (auto i = 0; i < values.currentPoint.size(); ++i) {
                        constrained[i] =
                                (values.currentPoint[i] == 0.0
                                 and values.currentGradient[i] > 0);
                        active = active || constrained[i];
                    }

                    constrained[2] = (std::abs(values.currentPoint[0] + values.currentPoint[1] - 1) <= eps)
                                     and (values.currentGradient[0] + values.currentGradient[1]) <= 0.0;
                    active = active || constrained[2];
                    if (not active)
                        return -values.currentHessian.inverse() * values.currentGradient;
                    if (
                            (constrained[0] and constrained[1])
                            or
                            (values.currentPoint[0] == 1. and values.currentPoint[1] == 0. and
                             (values.currentGradient[1] - values.currentGradient[0]) >= 0. and
                             values.currentGradient[0] <= 0.)
                            or
                            (values.currentPoint[0] == 0. and values.currentPoint[1] == 1. and
                             (values.currentGradient[1] - values.currentGradient[0]) <= 0. and
                             values.currentGradient[1] <= 0.)
                            ) {
                        // gradient is pushing outside the constrained area
                        return ZeroVec;
                    }

                    Matrix Hi = values.currentHessian.inverse();
                    if (constrained[0]) {
                        Hi.row(0).fill(0.);
                        Hi.col(0).fill(0.);
                    } else if (constrained[1]) {
                        Hi.row(1).fill(0.);
                        Hi.col(1).fill(0.);
                    } else if (constrained[2]) {
                        Hi = P3 * Hi * P3;
                    }
                    return -Hi * values.currentGradient;
                } else {
                    bool active = ((values.currentPoint[0] == 0.0
                                    and values.currentGradient[0] > 0) or
                                   (values.currentPoint[0] == 1.0
                                    and values.currentGradient[0] < 0));
                    if (active)
                        return ZeroVec;
                    else
                        return -values.currentHessian.inverse() * values.currentGradient;
                }
            }

            /*!
             * @brief The class is clonable
             *
             * @return A clone of myself wrapped into a unique pointer
             */
            virtual
            std::unique_ptr<DescentDirectionBase<PRDIM>>
            clone() const override { return std::make_unique<NewtonDirection>(*this); }

        private:
            static inline const Vector ZeroVec{Vector::Zero()};
            static inline const Matrix P3{Matrix::Identity() - 0.5 * Matrix::Ones()};
            static constexpr double eps = 100. * std::numeric_limits<double>::epsilon();
        };

/*!
 *  Implements the classic BFGS quasi-Newton algorithm.
 */
        class BFGSDirection : public DescentDirectionBase<PRDIM> {
        public:
            Vector operator()(const OptimizationCurrentValues &values) {
                // First time is a gradient step
                if (firstTime) {
                    auto const n = values.currentPoint.size();
                    H = Eigen::MatrixXd::Identity(n, n);
                    firstTime = false;
                    this->previousValues = values;
                    return -values.currentGradient;
                }

                Vector const &g = values.currentGradient;
                Vector yk = g - this->previousValues.currentGradient;
                Vector sk =
                        values.currentPoint - this->previousValues.currentPoint;
                auto const yks = yk.dot(sk);
                // Correct approximate Hessian only if we maintain sdp property if not keep
                // the old one
                if (yks > this->smallNumber * sk.norm() * yk.norm()) {
                    Vector Hs;
                    Hs = H * sk;
                    H += (yk * yk.transpose()) / yks - (Hs * Hs.transpose()) / (sk.dot(Hs));
                }
                this->previousValues = values;
                Vector d = H.fullPivLu().solve(-g);
                return d;
            }

            /*!
             * You need to reset if you run another problem or you start from a different
             * initial point.
             * @note This is done inside LineSearchSolver
             */
            void reset() {
                this->firstTime = true;
            };

            /*!
              * @brief The class is clonable
              *
              * @return A clone of myself wrapped into a unique pointer
              */
            virtual
            std::unique_ptr<DescentDirectionBase<PRDIM>>
            clone() const override { return std::make_unique<BFGSDirection>(*this); }


        private:
            OptimizationCurrentValues previousValues;
            Eigen::MatrixXd H;
            bool firstTime{true};
            double const smallNumber = std::sqrt(std::numeric_limits<double>::epsilon());
        };

/*!
 * Implements BFGS with the direct computation of the approximate inverse of
 * the Hessian.
 */
        class BFGSIDirection : public DescentDirectionBase<PRDIM> {
        public:
            Vector operator()(const OptimizationCurrentValues &values) {
                // First time is a gradient step
                if (firstTime) {
                    auto n = values.currentPoint.size();
                    H = Eigen::MatrixXd::Identity(n, n);
                    firstTime = false;
                    this->previousValues = values;
                    return -values.currentGradient;
                }

                Vector const &g = values.currentGradient;
                Vector yk = g - this->previousValues.currentGradient;
                Vector sk =
                        values.currentPoint - this->previousValues.currentPoint;
                auto const yks = yk.dot(sk);
                // Correct approximate Hessian only if we maintain sdp property if not keep
                // the old one
                if (yks > this->smallNumber * sk.norm() * yk.norm()) {
                    H += sk * sk.transpose() * (yks + yk.transpose() * H * yk) / (yks * yks) -
                         (H * yk * sk.transpose() + sk * yk.transpose() * H) / yks;
                }
                this->previousValues = values;
                Vector d = -H * g;
                return d;
            };

            void reset() {
                this->firstTime = true;
            };

        private:
            OptimizationCurrentValues previousValues;
            Eigen::MatrixXd H;
            bool firstTime{true};
            double const smallNumber = std::sqrt(std::numeric_limits<double>::epsilon());

            /*!
              * @brief The class is clonable
              *
              * @return A clone of myself wrapped into a unique pointer
              */
            virtual
            std::unique_ptr<DescentDirectionBase<PRDIM>>
            clone() const override { return std::make_unique<BFGSIDirection>(*this); }
        };

/*!
 * Bazrzilain-Borwein
 */
        class BBDirection : public DescentDirectionBase<PRDIM> {
        public:
            Vector operator()(OptimizationCurrentValues const &values) {
                // First time is a gradient step
                if (firstTime) {
                    firstTime = false;
                    this->previousValues = values;
                    return -values.currentGradient;
                }

                Vector const &g = values.currentGradient;
                Vector yk = g - this->previousValues.currentGradient;
                Vector sk =
                        values.currentPoint - this->previousValues.currentPoint;
                auto yks = yk.dot(sk);
                auto ykk = yk.dot(yk);
                this->previousValues = values;
                // Correct approximate Hessian only if we maintain sdp property if not keep
                // the old one
                if (yks > this->smallNumber * sk.norm() * yk.norm() and ykk > smallNumber) {
                    // I use a mix of the two possible strategies
                    return -0.5 * (yks / ykk + sk.dot(sk) / yks) * g;
                    // Strategy 1
                    // return -(yks/ykk)*g;
                    // Strategy 2
                    // return -(sk.dot(sk)/yks)*g;
                } else {
                    return -g;
                }
            };

            void
            reset() override {
                firstTime = true;
            };

            /*!
              * @brief The class is clonable
              *
              * @return A clone of myself wrapped into a unique pointer
              */
            virtual
            std::unique_ptr<DescentDirectionBase<PRDIM>>
            clone() const override { return std::make_unique<BBDirection>(*this); }


        private:
            OptimizationCurrentValues previousValues;
            bool firstTime{true};
            double const smallNumber = std::sqrt(std::numeric_limits<double>::epsilon());
        };

/*!
 * Non-linear CG with Polak-Ribiere formula
 */
        class CGDirection : public DescentDirectionBase<PRDIM> {
        public:
            Vector operator()(const OptimizationCurrentValues &values) {
                // First time is a gradient step
                if (firstTime) {
                    firstTime = false;
                    this->prevDk = -values.currentGradient;
                } else {
                    Vector const &gk1 = values.currentGradient;
                    Vector const &gk =
                            this->previousValues.currentGradient;
                    Vector const &dk = this->prevDk;
                    // Polak Ribiere formula
                    this->prevDk =
                            -gk1 + (gk1.dot(gk1 - gk) / (gk.dot(gk))) * dk; // store for next time
                    if (prevDk.dot(gk1) > 0)
                        prevDk =
                                -gk1; // check if direction is a descent if not go back to gradient
                }
                this->previousValues = values; // store for next time
                return this->prevDk;
            };

            void
            reset() override {
                firstTime = true;
            };

            /*!
              * @brief The class is clonable
              *
              * @return A clone of myself wrapped into a unique pointer
              */
            virtual
            std::unique_ptr<DescentDirectionBase<PRDIM>>
            clone() const override { return std::make_unique<CGDirection>(*this); }


        private:
            OptimizationCurrentValues previousValues;
            bool firstTime{true};
            //! I need to keep track of previous descent direction
            Vector prevDk;
        };
    };
} // namespace apsc

#endif /* EXAMPLES_SRC_LINESEARCH_DESCENTDIRECTIONS_HPP_ */

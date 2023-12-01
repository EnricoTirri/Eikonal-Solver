/*
 * DecenntDirectionFactory.hpp
 *
 *  Created on: Mar 24, 2022
 *      Author: forma
 */

#ifndef EXAMPLES_SRC_LINESEARCH_DESCENTDIRECTIONFACTORY_HPP_
#define EXAMPLES_SRC_LINESEARCH_DESCENTDIRECTIONFACTORY_HPP_

#include "Factory.hpp"
#include "DescentDirections.hpp"

namespace apsc {
    /*!
     * @brief Load descent directions in the factory.
     *
     * Look at the source code for details. The definition is in the source file.
     * Here, for simplicity, I do not use the constructor attribute to have the automatic registration.
     * The function returns the reference to the only factory in the code.
     */
    template<size_t PRDIM>
    struct DirectionFactory {
        /*!
        * The Factory of descent directions. an instance of my generic Factory
        */
        using DescentDirectionFactory = GenericFactory::Factory<DescentDirectionBase<PRDIM>, std::string>;

        DescentDirectionFactory &loadDirections() {
            // get the factory
            DescentDirectionFactory &theFactory = DescentDirectionFactory::Instance();
            //
            theFactory.add("GradientDirection",
                           []() { return std::make_unique<DescentDirection<PRDIM>::GradientDirection>(); });
            theFactory.add("BFGSDirection",
                           []() { return std::make_unique<DescentDirection<PRDIM>::BFGSDirection>(); });
            theFactory.add("BFGSIDirection",
                           []() { return std::make_unique<DescentDirection<PRDIM>::BFGSIDirection>(); });
            theFactory.add("BBDirection",
                           []() { return std::make_unique<DescentDirection<PRDIM>::BBDirection>(); });
            theFactory.add("CGDirection",
                           []() { return std::make_unique<DescentDirection<PRDIM>::CGDirection>(); });
            theFactory.add("NewtonDirection",
                           []() { return std::make_unique<DescentDirection<PRDIM>::NewtonDirection>(); });
            return theFactory;
        }
        /*!
         * Anothe possibility is to have the factory as global variable
         */
        //extern DescentDirectionFactory & descentDirectionFactory;
    };
}


#endif /* EXAMPLES_SRC_LINESEARCH_DESCENTDIRECTIONFACTORY_HPP_ */

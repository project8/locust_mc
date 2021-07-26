/*
 * LMCCavityModes.hh
 *
 *  Created on: Feb 28, 2020
 *      Author: pslocum
 */

#ifndef LMCCAVITYMODES_HH_
#define LMCCAVITYMODES_HH_

#include "LMCPowerCombiner.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class CavityModes
 @author P. Slocum
 @brief Derived class describing resonant cavity
 @details
 Available configuration options:
 No input parameters
 */


    class CavityModes: public PowerCombiner
    {

        public:
            CavityModes();
            virtual ~CavityModes();
            virtual bool Configure( const scarab::param_node& aNode );


        private:



    };


} /* namespace locust */

#endif /* LMCCAVITYMODES_HH_ */

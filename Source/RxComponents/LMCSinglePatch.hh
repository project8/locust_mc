/*
 * LMCSinglePatch.hh
 *
 *  Created on: Feb 28, 2020
 *      Author: pslocum
 */

#ifndef LMCSINGLEPATCH_HH_
#define LMCSINGLEPATCH_HH_

#include "LMCPowerCombinerParent.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class SinglePatch
 @author P. Slocum
 @brief Derived class describing a single patch.
 @details
 Available configuration options:
 No input parameters
 */


    class SinglePatch: public PowerCombinerParent
    {

        public:
            SinglePatch();
            virtual ~SinglePatch();
            virtual bool Configure( const scarab::param_node& aNode );

        	virtual bool SetVoltageDampingFactors();
        	virtual bool IsSinglePatch();



    };


} /* namespace locust */

#endif /* LMCSLOTTEDWAVEGUIDE_HH_ */

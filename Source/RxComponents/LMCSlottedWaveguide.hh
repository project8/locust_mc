/*
 * LMCSlottedWaveguide.hh
 *
 *  Created on: Feb 28, 2020
 *      Author: pslocum
 */

#ifndef LMCSLOTTEDWAVEGUIDE_HH_
#define LMCSLOTTEDWAVEGUIDE_HH_

#include "LMCPowerCombinerParent.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class SlottedWaveguide
 @author P. Slocum
 @brief Derived class describing the slot antenna
 @details
 Available configuration options:
 No input parameters
 */


    class SlottedWaveguide: public PowerCombinerParent
    {

        public:
            SlottedWaveguide();
            virtual ~SlottedWaveguide();
            virtual bool Configure( const scarab::param_node& aNode );

        	virtual bool SetVoltageDampingFactors();

        private:
        	double fImpedanceTransformation;  // != 1.0 only if not included in the HFSS TF.



    };


} /* namespace locust */

#endif /* LMCSLOTTEDWAVEGUIDE_HH_ */

/*
 * LMCWaveguideModes.hh
 *
 *  Created on: Mar 17, 2022
 *      Author: pslocum
 */

#ifndef LMCWAVEGUIDEMODES_HH_
#define LMCWAVEGUIDEMODES_HH_

#include "LMCPowerCombiner.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class WaveguideModes
 @author P. Slocum
 @brief Derived class describing summation at a probe of propagating waveguide modes excited by an electron
 inside the waveguide.
 @details
 Available configuration options:
 No input parameters
 */


    class WaveguideModes: public PowerCombiner
    {

        public:
            WaveguideModes();
            virtual ~WaveguideModes();
            virtual bool Configure( const scarab::param_node& aNode );
        	virtual bool AddOneModeToCavityProbe(Signal* aSignal, double excitationAmplitude, double dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex);


        private:
            double fVoltagePhase;



    };


} /* namespace locust */

#endif /* LMCWAVEGUIDEMODES_HH_ */

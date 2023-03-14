/*
 * LMCCavityModes.hh
 *
 *  Created on: Mar. 18, 2022
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
        	virtual bool AddOneModeToCavityProbe(Signal* aSignal, std::vector<double> particleXP, std::vector<double> excitationAmplitude, std::vector<double> EFieldAtProbe, std::vector<double> dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, bool initParticle);
        	virtual bool AddOneSampleToRollingAvg(int l, int m, int n, std::vector<double> excitationAmplitude, unsigned sampleIndex);
            double GetCavityProbeGain();
            void SetCavityProbeGain( double aGain );
            double GetVoltagePhase();
            void SetVoltagePhase( double aPhase );



        private:
            double fProbeGain;
            double fOrbitPhase;
            double fVoltagePhase;
            std::vector<std::vector<std::vector<double>>> fRollingAvg;
            std::vector<std::vector<std::vector<int>>> fCounter;




    };


} /* namespace locust */

#endif /* LMCCAVITYMODES_HH_ */

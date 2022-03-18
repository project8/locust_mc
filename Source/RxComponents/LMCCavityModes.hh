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
        	virtual bool AddOneModeToCavityProbe(Signal* aSignal, double excitationAmplitude, double dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex);
        	virtual bool AddOneSampleToRollingAvg(int l, int m, int n, double excitationAmplitude, unsigned sampleIndex);
            std::vector<double> GetCavityProbeZ();
            void SetCavityProbeZ ( std::vector<double> aVector );
            std::vector<double> GetCavityProbeTheta();
            void SetCavityProbeTheta ( std::vector<double> aVector );
            int GetNCavityProbes();
            void SetNCavityProbes( int aNumberOfProbes );
            double GetCavityProbeInductance();
            void SetCavityProbeInductance( double anInductance );
            bool SetCavityProbeLocations(int nCavityProbes, double cavityLength);


        private:
            int fNCavityProbes;
            double fCavityProbeInductance;
            std::vector<double> fCavityProbeZ;
            std::vector<double> fCavityProbeTheta;
            std::vector<std::vector<std::vector<double>>> fRollingAvg;
            std::vector<std::vector<std::vector<int>>> fCounter;




    };


} /* namespace locust */

#endif /* LMCCAVITYMODES_HH_ */

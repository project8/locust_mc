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
#ifdef ROOT_FOUND
    #include "LMCRootHistoWriter.hh"
#endif


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
        	virtual bool AddOneModeToCavityProbe(int l, int m, int n, Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, int channelIndex, bool initParticle);
        	virtual bool AddOneSampleToRollingAvg(int l, int m, int n, double excitationAmplitude, unsigned sampleIndex);
            std::vector<std::vector<std::vector<std::vector<double>>>> GetVoltagePhase(int aChannel, int l, int m, int n);
            void SetVoltagePhase( double aPhase, int aChannel, int l, int m, int n);
        	bool WriteRootHisto();



        private:
            double fOrbitPhase;
            std::vector<std::vector<std::vector<std::vector<double>>>> fVoltagePhase;
            std::vector<std::vector<std::vector<double>>> fRollingAvg;
            std::vector<std::vector<std::vector<int>>> fCounter;
            FILE *fp;
#ifdef ROOT_FOUND
            FileWriter* fRootHistoWriter;
#endif




    };


} /* namespace locust */

#endif /* LMCCAVITYMODES_HH_ */

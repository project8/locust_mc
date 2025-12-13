/*
 * LMCCavityModes.hh
 *
 *  Created on: Mar. 18, 2022
 *      Author: pslocum
 */

#ifndef LMCCAVITYMODES_HH_
#define LMCCAVITYMODES_HH_

#include "LMCKassLocustInterface.hh"
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
       - "norm-check":  bool [false] Record energy depositions of normalized modes into text
       file and Root histogram file.  Each record in the text/Root file represents an updated
       rolling average of relative mode energy depositions.

 */


    class CavityModes: public PowerCombiner
    {

        public:
            CavityModes();
            virtual ~CavityModes();
            virtual bool Configure( const scarab::param_node& aNode );
            virtual bool AddOneModeToCavityProbe(int mu, Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, int channelIndex, bool initParticle);
            double GetVoltagePhase(int aChannel, int mu);
            void SetVoltagePhase( double aPhase, int aChannel, int mu);
            bool SetChannelPhaseShifts();
            virtual bool SizeNChannels(int aNumberOfChannels);



        private:
            std::vector<std::vector<double>> fVoltagePhase;
            bool fModeMap;
            std::vector<std::vector<int>> fModeSet;
            std::vector<std::vector<double>> fChannelPhaseShifts;
            kl_interface_ptr_t fInterface;




    };


} /* namespace locust */

#endif /* LMCCAVITYMODES_HH_ */

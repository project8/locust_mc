/*
 * LMCCavitySignalGenerator.hh
 *
 *  Created on: Mar 30, 2021
 *      Author: pslocum
 */

#ifndef LMCCAVITYSIGNALGENERATOR_HH_
#define LMCCAVITYSIGNALGENERATOR_HH_

#include "LMCThreeVector.hh"
#include <boost/math/special_functions/bessel.hpp>
#include "LMCGenerator.hh"
#include "LMCCavityModes.hh" // : LMCPowerCombiner
#include "LMCKassLocustInterface.hh"
#include "LMCKassCurrentTransmitter.hh"
#include "LMCField.hh"
#include "LMCCylindricalCavity.hh" // : LMCField
#include "LMCFieldBuffer.hh"
#include <vector>
#include <sstream>
#include <string>
#include "LMCException.hh"
#include <algorithm>    // std::min



namespace locust
{

    /*!
     @class CavitySignalGenerator
     @author P. L. Slocum

     @brief Generate signal in cavity and detect it with antennas.

     @details
     Operates in [tbd] domain

     Configuration name: "cavity-signal"

     Available configuration options:
     - "param-name": type -- Description
     - "lo-frequency" : double -- local oscillator frequency
     - "xml-filename" : std::string -- the name of the xml locust config file.
     - "lo-frequency":  local oscillator frequency in Hz.

    */

    class CavitySignalGenerator : public Generator
    {
        public:

            CavitySignalGenerator( const std::string& aName = "cavity-signal" );
            virtual ~CavitySignalGenerator();

            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;
              
            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );
            std::vector<std::vector<std::vector<double>>> CalculateNormFactors(int nModes, bool TE);
            std::vector<int> ModeFilter(unsigned whichMode);
            void CheckNormalization();
            void PrintModeMaps();




        private:
            double fLO_Frequency;
            int fNModes;
            int fNPreEventSamples;  // spacing between events.  constant for now, could be randomized.
            int fThreadCheckTime;  // time (ms) to check for response from Kass thread.
            double fArrayRadius;
            std::string gxml_filename;
            bool fKassNeverStarted;
            bool fSkippedSamples;
            double fphiLO; // voltage phase of LO in radians;
                        std::vector<std::deque<double>> fStupidBuffer;
//            std::vector<std::deque<double>> fLocalFIRfrequencyBuffer;
//            std::vector<std::deque<double>> fLocalElementFIRBuffer;



            void ReadBesselZeroes(std::string filename, bool prime);
            void KassiopeiaInit(const std::string &aFile);
            void WakeBeforeEvent();
            bool ReceivedKassReady();
            bool DriveMode(Signal* aSignal, int nFilterBinsRequired, double dtFilter, unsigned index);
            double GetModeScalingFactor(std::vector<double> tKassParticleXP, int channelIndex);
            void InitializeBuffers(unsigned filterbuffersize);
            std::vector<double> GetCavityNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP);
            double GetCavityDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized);
            double GetCavityFIRSample(std::vector<double> tKassParticleXP, std::vector<std::deque<double>> tLocalFIRfrequencyBuffer, std::vector<std::deque<double>> tLocalElementFIRBuffer,int nFilterBinsRequired, double dtFilter);
            std::vector<std::deque<double>> create_copy(std::vector<std::deque<double>> const &vec);


            bool DoGenerate( Signal* aSignal );
            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );
            bool (CavitySignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

            PowerCombiner* fPowerCombiner;

            kl_interface_ptr_t fInterface;
            FILE *fp;

    };

} /* namespace locust */

#endif /* LMCCAVITYSIGNALGENERATOR_HH_ */

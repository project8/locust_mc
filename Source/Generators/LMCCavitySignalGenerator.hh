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
#include "LMCEquivalentCircuit.hh"
#include "LMCKassLocustInterface.hh"
#include "LMCKassCurrentTransmitter.hh"
#include "LMCField.hh"
#include "LMCCylindricalCavity.hh" // : LMCField
#include "LMCRectangularWaveguide.hh" // : LMCField
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
     - "cavity-radius" : double -- ideal cylindrical cavity radius.
     - "cavity-length" : double -- ideal cylindrical cavity length.
     - "n-modes" : int -- range of l, m, and n indices used to configure mode normalizations.
     	 	 However, presently the only mode being simulated is 011.
     - "lo-frequency" : double -- local oscillator frequency
     - "xml-filename" : std::string -- the name of the xml locust config file.
     - "lo-frequency":  local oscillator frequency in Hz.
     - "bypass-tf":  bool(false) -- if true, set FIR convolution output to 1.0
     - "norm-check": bool(false) -- if true, calculate weighted running averages of J \cdot E
      	 for all modes with indices of order < fNModes, and write the avgs to an intermediate file.
     - "equivalentR": Resistance from equivalent RLC circuit in Ohms
     - "equivalentL": Inductance from equivalent RLC circuit in Henries
     - "equivalentC": Capacitance from equivalent RLC circuit in Farads

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
            bool ModeSelect(int l, int m, int n, bool eGun);
            void CheckNormalization();
            void PrintModeMaps();
            double ScaleEPoyntingVector(double fcyc);
	    double ScaleAmplitude(int l, int m, int n, double fcyc);


        private:
            bool fE_Gun;
            double fLO_Frequency;
            int fNModes;
            int fNPreEventSamples;  // spacing between events.  constant for now, could be randomized.
            int fThreadCheckTime;  // time (ms) to check for response from Kass thread.
            double fArrayRadius;
            std::string gxml_filename;
            bool fKassNeverStarted;
            bool fSkippedSamples;
            double fphiLO; // voltage phase of LO in radians;
            bool fBypassTF;
            bool fNormCheck;



            void ReadBesselZeroes(std::string filename, bool prime);
            void KassiopeiaInit(const std::string &aFile);
            void WakeBeforeEvent();
            bool ReceivedKassReady();
            bool DriveMode(Signal* aSignal, int nFilterBinsRequired, double dtFilter, unsigned index);
            double GetModeScalingFactor(std::vector<double> tKassParticleXP, int channelIndex);
            void InitializeBuffers(unsigned filterbuffersize);
            std::vector<double> GetCavityNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP);
            std::vector<double> GetWaveguideNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP);
            double GetCavityDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized);
            double GetWaveguideDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized);
            double GetCavityFIRSample(std::vector<double> tKassParticleXP, std::vector<std::deque<double>> tLocalFIRfrequencyBuffer, std::vector<std::deque<double>> tLocalElementFIRBuffer,int nFilterBinsRequired, double dtFilter);


            bool DoGenerate( Signal* aSignal );
            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );
            bool (CavitySignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

            PowerCombiner* fPowerCombiner;
	    EquivalentCircuit* fEquivalentCircuit;

            kl_interface_ptr_t fInterface;
            FILE *fp;

    };

} /* namespace locust */

#endif /* LMCCAVITYSIGNALGENERATOR_HH_ */

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
#include "LMCWaveguideModes.hh" // : LMCPowerCombiner
#include "LMCEquivalentCircuit.hh" // : LMCAnalyticResponseFunction
#include "LMCDampedHarmonicOscillator.hh" // : LMCAnalyticResponseFunction
#include "LMCKassLocustInterface.hh"
#include "LMCKassCurrentTransmitter.hh"
#include "LMCFieldCalculator.hh"
#include "LMCField.hh"
#include "LMCCylindricalCavity.hh" // : LMCField
#include "LMCRectangularWaveguide.hh" // : LMCField
#include "LMCFIRFileHandler.hh"
#include "LMCTFFileHandler.hh"
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
     - "e-gun": Select e-gun configuration instead of cavity [false].
     - "center-to-short": distance [0.05 m] from center of e-gun waveguide to reflecting short.
     - "center-to-antenna": distance [0.05 m] from center of e-gun waveguide to antenna.
     - "waveguide-short":  optional presence/absence of reflecting short [true].
     - "back-reaction": optional waveguide back reaction in e-gun.  default to [true] if waveguide-short is present.
     - "direct-kass-power":  In e-gun, overrides calculated waveguide signal amplitudes and replaces
     	 them with sqrt(KassPower).  This is for cross-checking the more detailed signal calculations.
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



        private:
            double fLO_Frequency;
            double fDeltaT;
            int fNModes;
            int fNPreEventSamples;  // spacing between events.  constant for now, could be randomized.
            int fThreadCheckTime;  // time (ms) to check for response from Kass thread.
            std::string gxml_filename;
            bool fKassNeverStarted;
            bool fSkippedSamples;
            double fphiLO; // voltage phase of LO in radians;
            double fAvgDotProductFactor;
            double fdtFilter;
            bool fBypassTF;
            bool fNormCheck;
            bool fModeMaps;
            bool fTE; // (if false, use TM modes.)
            bool fIntermediateFile;
            bool fUseDirectKassPower;




            void ReadBesselZeroes(std::string filename, bool prime);
            void KassiopeiaInit(const std::string &aFile);
            void WakeBeforeEvent();
            bool ReceivedKassReady();
            bool DriveMode(Signal* aSignal, unsigned index);


            bool DoGenerate( Signal* aSignal );
            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );
            bool (CavitySignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

            PowerCombiner* fPowerCombiner;
            FieldCalculator* fFieldCalculator;
            TFReceiverHandler* fTFReceiverHandler;
            AnalyticResponseFunction* fAnalyticResponseFunction;

            kl_interface_ptr_t fInterface;
            FILE *fp;

    };

} /* namespace locust */

#endif /* LMCCAVITYSIGNALGENERATOR_HH_ */

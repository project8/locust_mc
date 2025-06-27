/*
 * LMCPowerNormCavitySignalGenerator.hh
 *
 *  Created on: Jun 27, 2025
 *      Author: rmueller
 */

#ifndef LMCPOWERNORMCAVITYSIGNALGENERATOR_HH_
#define LMCPOWERNORMCAVITYSIGNALGENERATOR_HH_

#include "LMCThreeVector.hh"
#include <boost/math/special_functions/bessel.hpp>
#include "LMCGenerator.hh"
#include "LMCCavityModes.hh" // : LMCPowerCombiner
#include "LMCWaveguideModes.hh" // : LMCPowerCombiner
#include "LMCEquivalentCircuit.hh" // : LMCAnalyticResponseFunction
#include "LMCDampedHarmonicOscillator.hh" // : LMCAnalyticResponseFunction
#include "LMCKassLocustInterface.hh"
#include "LMCKassCurrentTransmitter.hh"
#include "LMCPowerNormFieldCalculator.hh"
#include "LMCField.hh"
#include "LMCCylindricalCavity.hh" // : LMCField
#include "LMCRectangularCavity.hh" // : LMCField
#include "LMCRectangularWaveguide.hh" // : LMCField
#include "LMCFIRFileHandler.hh"
#include "LMCTFFileHandler.hh"
#include "LMCCavityUtility.hh"
#include "LMCAliasingUtility.hh"
#include "LMCDistributionInterface.hh"
#include <vector>
#include <sstream>
#include <string>
#include "LMCException.hh"
#include <stdlib.h>
#include <time.h>




namespace locust
{

    /*!
     @class PowerNormCavitySignalGenerator
     @author P. L. Slocum

     @brief Generate signal in cavity and detect it with antennas.

     @details
     Operates in time domain

     Configuration name: "cavity-signal"

     Available configuration options:
     - "xml-filename" : std::string -- the name of the xml locust config file.
     - "lo-frequency":  double -- local oscillator frequency in Hz.
     - "bypass-tf":  bool(false) -- if true, set FIR convolution output to 1.0 .
     - "back-reaction": bool [true] -- Default is true unless defined as false, or if "waveguide-short" = false.
     - "direct-kass-power":  bool [true] In waveguide simulation, populates signal amplitudes with sqrt(KassPower).
     This is useful for cross-checking the more detailed signal calculations.  Value changes to [false]
     if tf-receiver-filename is defined in the config file.
     - "tf-receiver-filename": string [""] Full path to HFSS output file.  An example is the file
     data/WEGA_Impedance_Center.txt.  If this file is specified, then the parameter tf-receiver-bin-width
     must also be specified.
     - "tf-receiver-bin-width" double [100.e6] Frequency bin width contained in tf-receiver-filename.
     TO-DO:  Calculate this automatically instead of having to specify this parameter.
    */

    class PowerNormCavitySignalGenerator : public Generator
    {
        public:

            PowerNormCavitySignalGenerator( const std::string& aName = "cavity-signal" );
            virtual ~PowerNormCavitySignalGenerator();

            bool Configure( const scarab::param_node& aNode );
            bool ConfigureInterface(Signal* aSignal);
            bool RecordRunParameters(Signal* aSignal);
            bool CrossCheckCavityConfig();
            bool CrossCheckAliasing(Signal* aSignal, double dopplerFrequency );

            void Accept( GeneratorVisitor* aVisitor ) const;
              
            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );
            void CheckNormalization();

            const scarab::param_node* GetParameters();
            void SetParameters( const scarab::param_node& aNode );
            bool SetSeed(int aSeed);
            int GetSeed();



        private:
            double fLO_Frequency;
            double fDeltaT;
            int fNPreEventSamples;  // spacing between events.
            int fTrackDelaySeed; // Seed to randomize fNPreEventSamples
            bool fRandomPreEventSamples;
            int fThreadCheckTime;  // time (ms) to check for response from Kass thread.
            std::string gxml_filename;
            bool fKassNeverStarted;
            bool fAliasedFrequencies;
            bool fOverrideAliasing;
            double fphiLO; // voltage phase of LO in radians;
            bool fBypassTF;
            bool fNormCheck;
            bool fUseDirectKassPower;
            bool fAliasingIsChecked;
            bool fUnitTestRootFile;
            std::vector<std::vector<int>> fModeSet;




            void KassiopeiaInit(const std::string &aFile);
            void WakeBeforeEvent();
            bool TryWakeAgain();
            bool ReceivedKassReady();
            bool DriveMode(Signal* aSignal, unsigned index);
            bool RandomizeStartDelay();


            bool DoGenerate( Signal* aSignal );
            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );
            bool (PowerNormCavitySignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

            PowerCombiner* fPowerCombiner;
            FieldCalculator* fFieldCalculator;
            TFReceiverHandler* fTFReceiverHandler;
            AnalyticResponseFunction* fAnalyticResponseFunction;

            kl_interface_ptr_t fInterface;
            FILE *fp;

            const scarab::param_node* fParam;

            std::shared_ptr< BaseDistribution> fStartDelayDistribution;
            DistributionInterface fDistributionInterface;


    };

} /* namespace locust */

#endif /* LMCPOWERNORMCAVITYSIGNALGENERATOR_HH_ */

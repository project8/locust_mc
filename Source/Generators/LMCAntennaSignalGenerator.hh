/*
 * LMCAntennaSignalGenerator.hh
 *
 *  Created on: May 4, 2019
 *      Author: pslocum
 */

#ifndef LMCANTENNASIGNALGENERATOR_HH_
#define LMCANTENNASIGNALGENERATOR_HH_

#include "LMCThreeVector.hh"
#include "LMCGenerator.hh"
#include "LMCChannel.hh"
#include "LMCPatchAntenna.hh"
#include "LMCPowerCombiner.hh"
#include "LMCFieldEstimator.hh"


namespace locust
{

    /*!
     @class AntennaSignalGenerator
     @author P. L. Slocum

     @brief Generate artibtrary plane wave for calibration and detect with patch array.

     @details
     Operates in time space

     Configuration name: "antenna-signal"

     Available configuration options:
     - "param-name": type -- Description
     - "lo-frequency" : double -- local oscillator frequency
     - "antenna-frequency": 25.9281e9,
     - "array-radius": 0.05,
     - "npatches-per-strip": 22,
     - "patch-spacing": 0.0068,
     - "feed": "one-quarter",  -- power-combining scheme.
     - "phase-delay": 1, -- consider phase delay between patches?
     - "voltage-damping": 1 -- consider power-combining scheme in "feed"?
     - "input-signal-type": string -- By default considers sin wave, could change based on requirements. 
          
    */
    class AntennaSignalGenerator : public Generator
    {
        public:

            AntennaSignalGenerator( const std::string& aName = "antenna-signal" );
            virtual ~AntennaSignalGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;
              
            /*void AddOnePatchVoltageToStripSum(Signal* aSignal, double VoltageAmplitude, double VoltagePhase, double phi_LO, unsigned channelindex, unsigned z_index, double DopplerFrequency);
            double GetMismatchFactor(double f);
            double GetAOIFactor(double AOI, LMCThreeVector PatchNormalVector);
            double GetVoltageAmpFromAntenna();
	    */

        private:
	    FieldEstimator fFieldEstimator;
	    int fInputSignalType;
	    double fInputFrequency;// in GHz
            double fInputAmplitude;
           /* std::vector< Channel<PatchAntenna> > allChannels; //Vector that contains pointer to all channels
            std::vector<LMCThreeVector > rReceiver; //Vector that contains 3D position of all points at which the fields are evaluated (ie. along receiver surface)
            double fLO_Frequency;  // typically defined by a parameter in json file.
            double fRF_Frequency;  // typically defined by a parameter in json file.
            double fRF_Amplitude;  // typically defined by a parameter in json file.
            double fArrayRadius;  // from json file.
            int fNPatchesPerStrip; // from json file.
            double fPatchSpacing; // from json file.
            int fPowerCombiner;
            bool fVoltageDamping;*/

            void GenerateSignal(Signal* );

            bool DoGenerate( Signal* aSignal );
            //void* DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal);
            //void InitializePatchArray();


            //double phiLO_t; // voltage phase of LO in radians;
            //double VoltagePhase_t[10000];

    };

} /* namespace locust */

#endif /* LMCANTENNASIGNALGENERATOR_HH_ */

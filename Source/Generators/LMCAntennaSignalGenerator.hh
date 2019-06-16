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
              
        private:
	    FieldEstimator fFieldEstimator;
	    int fInputSignalType;
	    double fInputFrequency;// in GHz
            double fInputAmplitude;

            void GenerateSignal(Signal* );

            bool DoGenerate( Signal* aSignal );

    };

} /* namespace locust */

#endif /* LMCANTENNASIGNALGENERATOR_HH_ */

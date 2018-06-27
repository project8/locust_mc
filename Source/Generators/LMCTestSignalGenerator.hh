/*
 * LMCTestSignalGenerator.hh
 *
 *  Created on: Jan 14 2015
 *      Author: plslocum and nsoblath
 */

#ifndef LMCTESTSIGNALGENERATOR_HH_
#define LMCTESTSIGNALGENERATOR_HH_

#include "LMCGenerator.hh"
#include "LMCRunLengthCalculator.hh"


namespace scarab
{
  class param_node;
}

namespace locust
{
  class Digitizer;

    /*!
     @class TestSignalGenerator
     @author P. L. Slocum

     @brief Add Sine Wave to the signal.

     @details
     Can operate in time or frequency space

     Configuration name: "test-signal"

     Available configuration options:
     - "frequency": double -- Frequency of the sine wave.
     - "amplitude": double -- Amplitude of the sine wave.
     - "domain": string -- Determines whether the sinusoidal test signal is generated in the time 
            or frequency domain
    
     Available options: "time" and "freq" [default]

    */
    class TestSignalGenerator : public Generator
    {
        public:
            TestSignalGenerator( const std::string& aName = "test-signal" );
            virtual ~TestSignalGenerator();

            bool Configure( const scarab::param_node* aNode );
      bool Configure2( const Digitizer* aDig );


            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetFrequency() const;
            void SetFrequency( double aFrequency );

            double GetAmplitude() const;
            void SetAmplitude( double aAmplitude );

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );


        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (TestSignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

            double fLO_frequency;
            double fAmplitude;

            
    };

} /* namespace locust */

#endif /* LMCTestSignalGENERATOR_HH_ */


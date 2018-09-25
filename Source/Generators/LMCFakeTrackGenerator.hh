/*
 * LMCFakeTrackGenerator.hh
 *
 *  Created on: Aug 8 2018
 *      Author: plslocum
 */

#ifndef LMCFAKETRACKGENERATOR_HH_
#define LMCFAKETRACKGENERATOR_HH_

#include "LMCGenerator.hh"
#include "LMCRunLengthCalculator.hh"
#include <random>


namespace scarab
{
  class param_node;
}

namespace locust
{
  class Digitizer;

    /*!
     @class FakeTrackGenerator
     @author P. L. Slocum

     @brief Add Sine Wave to the signal.

     @details
     Operates in time.

     Configuration name: "fake-track"

     Available configuration options:
     - "frequency": double -- Frequency of the sine wave.
     - "amplitude": double -- Amplitude of the sine wave.
     - "domain": string -- Determines whether the sinusoidal test signal is generated in the time 
            or frequency domain
    

    */
    class FakeTrackGenerator : public Generator
    {
        public:
            FakeTrackGenerator( const std::string& aName = "fake-track" );
            virtual ~FakeTrackGenerator();

            bool Configure( const scarab::param_node* aNode );
      bool Configure2( const Digitizer* aDig );


            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetFrequency() const;
            void SetFrequency( double aFrequency );

            double GetSignalPower() const;
            void SetSignalPower( double aPower );

            double GetStartFrequency() const;
            void SetStartFrequency( double aFrequency );

            double GetStartVPhase() const;
            void SetStartVPhase( double aPhase );

            double GetSlope() const;
            void SetSlope( double aSlope );

            double GetStartTime() const;
            void SetStartTime( double aTime );

            double GetEndTime() const;
            void SetEndTime( double aTime );


            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );


        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (FakeTrackGenerator::*fDoGenerateFunc)( Signal* aSignal );

            double fSignalPower;
            double fStartFrequency;
            double fStartVPhase;
            double fSlope;
            double fStartTime;
            double fEndTime;
            double fLO_frequency;

            

    };

} /* namespace locust */

#endif /* LMCFakeTrackGENERATOR_HH_ */


/*
 * LMCFakeTrackSignalGenerator.hh
 *
 *  Created on: Aug 8 2018
 *      Author: plslocum, L. Saldana
 */

#ifndef LMCFAKETRACKSIGNALGENERATOR_HH_
#define LMCFAKETRACKSIGNALGENERATOR_HH_

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
      @class FakeTrackSignalGenerator
      @author P. L. Slocum, L. Saldana

      @brief Generate custom, aka "fake", CRES track.

      @details
      Operates in time.

      Configuration name: "fake-track"

      Available configuration options:
      - "signal-power": double -- PSD of signal (W/Hz).
      - "start-frequency-max": double -- Upper bound for start frequency of signal (Hz); distribution: uniform.
      - "start-frequency-min": double -- Lower bound for start frequency of signal (Hz); distribution: uniform.
      - "track-length-mean": double -- Average of track length (s); distribution: exponential.
      - "start-vphase": double -- Starting voltage phase (V).
      - "slope-mean": double -- Mean value of Gaussian slope distribution (MHz/ms); distribution: gaussian.
      - "slope-std": double -- Standard deviation of Gaussian slope distribution (MHz/ms); distribution: gaussian.
      - "lo-frequency": double -- Frequency of local oscillator (Hz).
      - "start-time-max": double -- Upper bound for track start time (s); distribution: uniform.
      - "start-time-min": double -- Lower bound for track start time (s); distribution: uniform.
      - "ntracks-mean": double -- Average number of tracks per event (integer); distribution: exponential.
      - "random-seed": integer -- integer seed for random number generator for above pdfs, if set to 0 random_device will be used. 


*/
    class FakeTrackSignalGenerator : public Generator
    {
        public:
            FakeTrackSignalGenerator( const std::string& aName = "fake-track" );
            virtual ~FakeTrackSignalGenerator();

            bool Configure( const scarab::param_node* aNode );
            bool Configure2( const Digitizer* aDig );


            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetFrequency() const;
            void SetFrequency( double aFrequency );

            double GetSignalPower() const;
            void SetSignalPower( double aPower );

            double GetStartFrequencyMax() const;
            void SetStartFrequencyMax( double aFrequencyMax );

            double GetStartFrequencyMin() const;
            void SetStartFrequencyMin( double aFrequencyMin );

            double GetTrackLengthMean() const;
            void SetTrackLengthMean( double aTrackLengthMean );

            double GetStartVPhase() const;
            void SetStartVPhase( double aPhase );

            double GetSlopeMean() const;
            void SetSlopeMean( double aSlopeMean );

            double GetSlopeStd() const;
            void SetSlopeStd( double aSlopeStd );

            double GetStartTimeMax() const;
            void SetStartTimeMax( double aTimeMax );

            double GetStartTimeMin() const;
            void SetStartTimeMin( double aTimeMin );

            double GetNTracksMean() const;
            void SetNTracksMean( double aNTracksMean );

            int GetRandomSeed() const;
            void SetRandomSeed(  int aRandomSeed );


            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );


        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (FakeTrackSignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

            double fSignalPower;
            double fStartFrequencyMax;
            double fStartFrequencyMin;
            double fStartVPhase;
            double fSlopeMean;
            double fSlopeStd;
            double fStartTimeMax;
            double fStartTimeMin;
            double fLO_frequency;
            double fTrackLengthMean;
            double fNTracksMean;
            int fRandomSeed;


    };

} /* namespace locust */

#endif /* LMCFAKETRACKSIGNALGENERATOR_HH_ */


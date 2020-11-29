/*
 * LMCFakeFreeSpaceSignalGenerator.hh
 *
 *  Created on: Nov. 26, 2020
 *      Author: buzinsky
 */

#ifndef LMCFAKEFREESPACESIGNALGENERATOR_HH_
#define LMCFAKEFREESPACESIGNALGENERATOR_HH_

#include "LMCGenerator.hh"
#include "LMCEvent.hh"
#include "LMCRunParameters.hh"
#include "LMCRootTreeWriter.hh"
#include "LMCDistributionInterface.hh"

#include "LMCChannel.hh"
#include "LMCPatchAntenna.hh"


#include <random>
#include <vector>


namespace scarab
{
    class param_node;
}

namespace locust
{

    /*!
     @class FakeFreeSpaceSignalGenerator
     @author N. Buzinsky

     @brief Generate custom, aka "fake" Free-Space CRES events.

     @details
     Operates in time space

     Configuration name: "free-space-fake"

     Available configuration options:
      - "signal-power": double -- PSD of signal (at 90 degrees) (W/Hz).
      - "start-vphase": double -- Starting voltage phase (V).
      - "lo-frequency": double -- Frequency of local oscillator (Hz).
      - "start-time-max": double -- Upper bound for track start time (s); distribution: uniform.
      - "start-time-min": double -- Lower bound for track start time (s); distribution: uniform.
      - "magnetic-field": double -- Magnetic field used to convert from frequency to energy (for jumpsize) (T).
      - "n-events": int -- Number of events per simulation, spaced by 0.5 ms (hardcoded).
      - "random-seed": integer -- integer seed for random number generator for above pdfs, if set to 0 random_device will be used.
      - "root-filename": string -- Name of output Root file.  This can have the same name as other

    */

    class FakeFreeSpaceSignalGenerator : public Generator
    {
        public:
            FakeFreeSpaceSignalGenerator( const std::string& aName = "free-space-fake" );
            virtual ~FakeFreeSpaceSignalGenerator();

            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetLOFrequency() const;
            void SetLOFrequency( double aFrequency );

            double GetSignalPower() const;
            void SetSignalPower( double aPower );

            double GetStartVPhase() const;
            void SetStartVPhase( double aPhase );

            double GetStartTimeMax() const;
            void SetStartTimeMax( double aTimeMax );

            double GetStartTimeMin() const;
            void SetStartTimeMin( double aTimeMin );

            double GetNTracksMean() const;
            void SetNTracksMean( double aNTracksMean );

            double GetBField() const;
            void SetBField( double aBField );

            double GetAntennaRadius() const;
            void SetAntennaRadius( double aAntennaRadius );

            double GetGradBFrequency() const;
            void SetGradBFrequency( double aGradBFrequency );

            int GetRandomSeed() const;
            void SetRandomSeed(  int aRandomSeed );

            int GetNEvents() const;
            void SetNEvents(  int aNEvents );

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );

            double rel_cyc(const double& aEnergy, const double& aMagneticField) const;
            double rel_energy(const double& aFrequency, const double& aMagneticField) const;

            LMCThreeVector GetElectronPosition(const double& aRadialPhase);

            double GetCRESAmplitude(const LMCThreeVector& aElectronPosition, const unsigned channelIndex  );
            double GetCRESPhase(const double& aElectronTime, const LMCThreeVector& aElectronPosition, const unsigned& channelIndex  );

            void SetTrackProperties(Track &aTrack, int TrackID, double aTimeOffset);
            void InitiateEvent(Event* anEvent, int eventID);

            double fSlope;
            double fTrackLength;
            double fStartTime;
            double fEndTime;
            double fStartFrequency;
            double fCurrentFrequency;
            double fRadius;
            double fStartRadialPhase;
            double fCurrentRadialPhase;
            double fGradBFrequency;
            int fNTracks;


        private:
            bool DoGenerate( Signal* aSignal );
            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (FakeFreeSpaceSignalGenerator::*fDoGenerateFunc)( Signal* aSignal );

            void InitializeAntennaArray();

            std::vector< Channel<PatchAntenna> > allChannels; //Vector that contains pointer to all channels


            double fSignalPower;
            double fStartVPhase;

            std::shared_ptr< BaseDistribution> fStartEnergyDistribution;
            std::shared_ptr< BaseDistribution> fStartFrequencyDistribution;
            std::shared_ptr< BaseDistribution> fSlopeDistribution;
            std::shared_ptr< BaseDistribution> fRadiusDistribution;
            std::shared_ptr< BaseDistribution> fRadialPhaseDistribution;
            std::shared_ptr< BaseDistribution> fTrackLengthDistribution;

            double fStartTimeMax;
            double fStartTimeMin;
            double fAntennaRadius;
            double fLO_frequency;
            double fNTracksMean;
            double fBField;
            int fRandomSeed;
            int fNEvents;


            std::string fRootFilename;
            std::default_random_engine fRandomEngine;
            bool fUseEnergyDistribution;
            bool fUseFrequencyDistribution;

            DistributionInterface fDistributionInterface;
    };

} /* namespace locust */

#endif /* LMCFAKEFREESPACESIGNALGENERATOR_HH_ */

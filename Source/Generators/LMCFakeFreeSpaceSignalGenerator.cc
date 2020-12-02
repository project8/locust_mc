/*
 * LMCFakeFreeSpaceSignalGenerator.cc
 *
 *  Created on: Nov. 26, 2020
 *      Author: buzinsky
 */

#include "LMCFakeFreeSpaceSignalGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"
#include "LMCThreeVector.hh"
#include "path.hh"
#include <sstream>
#include <random>
#include <math.h>
#include <string>
#include <iostream>

using std::string;

namespace locust
{
    LOGGER( lmclog, "FakeFreeSpaceSignalGenerator" );

    MT_REGISTER_GENERATOR(FakeFreeSpaceSignalGenerator, "free-space-fake");

    FakeFreeSpaceSignalGenerator::FakeFreeSpaceSignalGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &FakeFreeSpaceSignalGenerator::DoGenerateTime ),
            fSignalPower( 0. ),
            fStartVPhase( 0. ),
            fStartTimeMin( 0. ),
            fStartTimeMax( 0. ),
            fLO_frequency( 0. ),
            fNTracksMean( 0. ),
            fBField( 1.0 ),
            fRandomSeed( 0 ),
            fNEvents( 1 ),
            fAntennaRadius( 0.05 ),
            fRadius( 0. ),
            fGradBFrequency( 0. ),
            fStartRadialPhase( 0. ),
            fCurrentRadialPhase( 0. ),
            fRandomEngine( 0 ),
            fRootFilename( "LocustEvent.root" ),
            fSlope( 0. ),
            fTrackLength( 0. ),
            fStartTime( 0. ),
            fEndTime( 0. ),
            fStartFrequency( 0. ),
            fCurrentFrequency( 0. ),
            fUseEnergyDistribution(false),
            fUseFrequencyDistribution(false),
            fNTracks(0)

    {
        fRequiredSignalState = Signal::kTime;
    }

    FakeFreeSpaceSignalGenerator::~FakeFreeSpaceSignalGenerator()
    {
    }

    bool FakeFreeSpaceSignalGenerator::Configure( const scarab::param_node& aParam )
    {

        if( aParam.has( "start-frequency" ) )
        {
            fStartFrequencyDistribution = fDistributionInterface.get_dist(aParam["start-frequency"].as_node());
            fUseFrequencyDistribution = true;
        }
        if( aParam.has( "start-energy" ) )
        {
            fStartEnergyDistribution = fDistributionInterface.get_dist(aParam["start-energy"].as_node());
            fUseEnergyDistribution = true;
        }

        if( aParam.has( "slope" ) )
        {
            fSlopeDistribution = fDistributionInterface.get_dist(aParam["slope"].as_node());
        }
        else
        {
            LWARN( lmclog, "Using default distribution: Slope = 0 ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            fSlopeDistribution = fDistributionInterface.get_dist(default_setting);
        }

        if( aParam.has( "radius" ) )
        {
            fRadiusDistribution = fDistributionInterface.get_dist(aParam["radius"].as_node());
        }
        else
        {
            LWARN( lmclog, "Using default distribution: radius = 0 ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            fRadiusDistribution = fDistributionInterface.get_dist(default_setting);
        }

        if( aParam.has( "radial-phase" ) )
        {
            fRadialPhaseDistribution = fDistributionInterface.get_dist(aParam["radial-phase"].as_node());
        }
        else
        {
            LWARN( lmclog, "Using default distribution: radial-phase = 0 ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            fRadialPhaseDistribution = fDistributionInterface.get_dist(default_setting);
        }

        if( aParam.has( "track-length" ) )
        {
            fTrackLengthDistribution = fDistributionInterface.get_dist(aParam["track-length"].as_node());
        }
        else
        {
            LWARN( lmclog, "Using default distribution: Track Length = 1e-4 ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            default_setting.add("value","1e-4");
            fTrackLengthDistribution = fDistributionInterface.get_dist(default_setting);
        }

        if( aParam.has( "signal-power" ) )
            SetSignalPower( aParam.get_value< double >( "signal-power", fSignalPower ) );

        if( aParam.has( "start-vphase" ) )
            SetStartVPhase( aParam.get_value< double >( "start-vphase", fStartVPhase ) );

        if( aParam.has( "start-time-max" ) )
            SetStartTimeMax( aParam.get_value< double >( "start-time-max", fStartTimeMax ) );

        if( aParam.has( "start-time-min" ) )
            SetStartTimeMin( aParam.get_value< double >( "start-time-min", fStartTimeMin ) );

        if( aParam.has( "lo-frequency" ) )
            SetLOFrequency( aParam.get_value< double >( "lo-frequency", fLO_frequency ) );

        if (aParam.has( "ntracks-mean") )
            SetNTracksMean( aParam.get_value< double >( "ntracks-mean",fNTracksMean) );

        if (aParam.has("magnetic-field") )
            SetBField(  aParam.get_value< double >("magnetic-field", fBField) );

        if( aParam.has( "antenna-radius" ) )
            SetAntennaRadius( aParam.get_value< double >( "antenna-radius", fAntennaRadius ) );

        if( aParam.has( "grad-B-frequency" ) )
            SetGradBFrequency( aParam.get_value< double >( "grad-B-frequency", fGradBFrequency ) );
        if (aParam.has( "random-seed") )
            SetRandomSeed(  aParam.get_value< int >( "random-seed",fRandomSeed) );

        if (aParam.has( "n-events") )
            SetNEvents(  aParam.get_value< int >( "n-events",fNEvents) );

        if( aParam.has( "root-filename" ) )
        {
            fRootFilename = aParam["root-filename"]().as_string();
        }

        if(fRandomSeed)
            fRandomEngine.seed(fRandomSeed);
        else
        {
            std::random_device rd;
            fRandomEngine.seed(rd());
        }

        fDistributionInterface.SetSeed(fRandomSeed);

        if( fUseFrequencyDistribution && fUseEnergyDistribution)
            LERROR( lmclog, "User specified both start frequency and start energy distribution! Please specify only one!");

        if( ! (fUseFrequencyDistribution || fUseEnergyDistribution))
        {
            LWARN( lmclog, "Using default distribution: Frequency = fLO + 50 MHz ");
            scarab::param_node default_setting;
            default_setting.add("name","dirac");
            default_setting.add("value",fLO_frequency + 50e6);
            fStartFrequencyDistribution = fDistributionInterface.get_dist(default_setting);
            fUseFrequencyDistribution = true;
        }

        if( aParam.has( "domain" ) )
        {
            string domain = aParam["domain"]().as_string();
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
                return false;
            }
        }


        return true;
    }

    void FakeFreeSpaceSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double FakeFreeSpaceSignalGenerator::GetSignalPower() const
    {
        return fSignalPower;
    }

    void FakeFreeSpaceSignalGenerator::SetSignalPower( double aPower )
    {
        fSignalPower = aPower;
        return;
    }

    double FakeFreeSpaceSignalGenerator::GetStartVPhase() const
    {
        return fStartVPhase;
    }

    void FakeFreeSpaceSignalGenerator::SetStartVPhase( double aPhase )
    {
        fStartVPhase = aPhase;
        return;
    }

    double FakeFreeSpaceSignalGenerator::GetStartTimeMin() const
    {
        return fStartTimeMin;
    }

    void FakeFreeSpaceSignalGenerator::SetStartTimeMin( double aTimeMin )
    {
        fStartTimeMin = aTimeMin;
        return;
    }

    double FakeFreeSpaceSignalGenerator::GetStartTimeMax() const
    {
        return fStartTimeMax;
    }

    void FakeFreeSpaceSignalGenerator::SetStartTimeMax( double aTimeMax )
    {
        fStartTimeMax = aTimeMax;
        return;
    }

    double FakeFreeSpaceSignalGenerator::GetLOFrequency() const
    {
        return fLO_frequency;
    }

    void FakeFreeSpaceSignalGenerator::SetLOFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
        return;
    }

    double FakeFreeSpaceSignalGenerator::GetNTracksMean() const
    {
        return fNTracksMean;
    }

    void FakeFreeSpaceSignalGenerator::SetNTracksMean( double aNTracksMean )
    {
        fNTracksMean = aNTracksMean;
        return;
    }


    double FakeFreeSpaceSignalGenerator::GetBField() const
    {
        return fBField;
    }

    void FakeFreeSpaceSignalGenerator::SetBField( double aBField )
    {
        fBField = aBField;
        return;
    }

    double FakeFreeSpaceSignalGenerator::GetAntennaRadius() const
    {
        return fAntennaRadius;
    }

    void FakeFreeSpaceSignalGenerator::SetAntennaRadius( double aAntennaRadius )
    {
        fAntennaRadius = aAntennaRadius;
        return;
    }

    double FakeFreeSpaceSignalGenerator::GetGradBFrequency() const
    {
        return fGradBFrequency;
    }

    void FakeFreeSpaceSignalGenerator::SetGradBFrequency( double aGradBFrequency )
    {
        fGradBFrequency = aGradBFrequency;
        return;
    }

    int FakeFreeSpaceSignalGenerator::GetRandomSeed() const
    {
        return fRandomSeed;
    }

    void FakeFreeSpaceSignalGenerator::SetRandomSeed( int aRandomSeed )
    {
        fRandomSeed = aRandomSeed;
        return;
    }
    
    int FakeFreeSpaceSignalGenerator::GetNEvents() const
    {
        return fNEvents;
    }

    void FakeFreeSpaceSignalGenerator::SetNEvents( int aNEvents )
    {
        fNEvents = aNEvents;
        return;
    }

    Signal::State FakeFreeSpaceSignalGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void FakeFreeSpaceSignalGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &FakeFreeSpaceSignalGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &FakeFreeSpaceSignalGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }

    double FakeFreeSpaceSignalGenerator::rel_cyc(const double& aEnergy, const double& aMagneticField) const
    {
        double tCycFreq = LMCConst::Q() * aMagneticField / LMCConst::M_el_kg();
        return tCycFreq / ( 1. + (aEnergy/LMCConst::M_el_eV()) ) / (2. * LMCConst::Pi()); // takes energy in eV, magnetic field in T, returns in Hz
    }

    double FakeFreeSpaceSignalGenerator::rel_energy(const double& aFrequency, const double& aMagneticField) const
    {
        double tGamma = LMCConst::Q() * aMagneticField / (2. * LMCConst::Pi() * aFrequency * LMCConst::M_el_kg()); //omega = q B / m Gamma
        return (tGamma - 1.) *LMCConst::M_el_eV(); // takes frequency in Hz, magnetic field in T, returns in kinetic energy eV
    }

    LMCThreeVector FakeFreeSpaceSignalGenerator::GetElectronPosition(const double& aRadialPhase)
    {
        return LMCThreeVector({fRadius * cos(aRadialPhase), fRadius * sin(aRadialPhase),0.});
    }



    double FakeFreeSpaceSignalGenerator::GetCRESPhase(const double& aElectronTime, const LMCThreeVector& aElectronPosition, const unsigned& channelIndex  )
    {
        PatchAntenna *currentPatch;
        currentPatch = &allChannels[channelIndex][0];  // only 1 slot / patch per channel
        
        LMCThreeVector tElectronAntenna = currentPatch->GetPosition() - aElectronPosition;

        double rRadius = tElectronAntenna.Magnitude();
        double tCurrentTheta = tElectronAntenna.AzimuthalAngle();

        double tCurrentAngularFrequency = 2. * LMCConst::Pi() * fCurrentFrequency;

        double tRetarded = aElectronTime - rRadius / LMCConst::C();
        return tCurrentAngularFrequency * tRetarded  - tCurrentTheta + fStartVPhase;
        

    }


    double FakeFreeSpaceSignalGenerator::GetCRESAmplitude(const LMCThreeVector& aElectronPosition, const unsigned channelIndex  ) 
    {
        PatchAntenna *currentPatch;
        currentPatch = &allChannels[channelIndex][0];  // only 1 slot / patch per channel
        const double tResistance = 50.; //Default electronic resistance antenna array (ohms)
        
        LMCThreeVector tElectronAntenna = currentPatch->GetPosition() - aElectronPosition;
        double rRadius = tElectronAntenna.Magnitude();

        return sqrt(tResistance * fSignalPower) * fAntennaRadius  / rRadius;
    }

    void FakeFreeSpaceSignalGenerator::InitializeAntennaArray()
    {
        const unsigned nChannels = fNChannels;
        const double patchRadius = fAntennaRadius;
        double theta;
        const double dThetaArray = 2. * LMCConst::Pi() / nChannels; //Divide the circle into nChannels

        PatchAntenna modelPatch;

        allChannels.resize(nChannels);

        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            theta = channelIndex * dThetaArray;

            modelPatch.SetCenterPosition({patchRadius * cos(theta) , patchRadius * sin(theta) , 0. }); 
            modelPatch.SetPolarizationDirection({sin(theta), -cos(theta), 0.}); 
            modelPatch.SetNormalDirection({-cos(theta), -sin(theta), 0.}); //Say normals point inwards
            allChannels[channelIndex].AddReceiver(modelPatch);
        }
    }

    void FakeFreeSpaceSignalGenerator::SetTrackProperties(Track &aTrack, int TrackID, double aTimeOffset)
    {
        double current_energy = 0.;
        double energy_loss = 0.;
        double new_energy = 0.;

        std::uniform_real_distribution<double> starttime_distribution(fStartTimeMin,fStartTimeMax);
        std::uniform_real_distribution<double> dist(0,1);

        if(TrackID==0)
        {
            if (aTimeOffset==0) // first event
        	{
                fStartTime = starttime_distribution(fRandomEngine);
        	}
            else
            {
                fStartTime = aTimeOffset;
            }

            if(fUseFrequencyDistribution)
            {
                fStartFrequency = fStartFrequencyDistribution->Generate();
            }
            else
            {
                fStartFrequency = rel_cyc(fStartEnergyDistribution->Generate(), fBField);
            }

            fRadius = fRadiusDistribution->Generate();
            fStartRadialPhase = fRadialPhaseDistribution->Generate();

            aTrack.Radius = fRadius;
            aTrack.RadialPhase = fStartRadialPhase;
            aTrack.StartTime = fStartTime;
            aTrack.StartFrequency = fStartFrequency;
        }

       else
       {
            fStartTime = fEndTime + 0.;  // old track endtime + margin=0.
            current_energy = rel_energy(fCurrentFrequency, fBField); // convert current frequency to energy

            energy_loss = 10.; //hardcoded jump (for now)
            new_energy = current_energy - energy_loss; // new energy after loss, in eV
            fStartFrequency = rel_cyc(new_energy, fBField);
            fCurrentFrequency = fStartFrequency;
            fStartRadialPhase = fCurrentRadialPhase;
            aTrack.StartTime = fEndTime + 0.; // margin of time is 0.
            aTrack.StartFrequency = fStartFrequency;
            aTrack.RadialPhase = fStartRadialPhase;
        }

        fSlope = fSlopeDistribution->Generate();

        fTrackLength = fTrackLengthDistribution->Generate();
        fEndTime = fStartTime + fTrackLength;  // reset endtime.
        aTrack.Slope = fSlope;
        aTrack.TrackLength = fTrackLength;
        aTrack.EndTime = aTrack.StartTime + aTrack.TrackLength;
        aTrack.LOFrequency = fLO_frequency;
        aTrack.TrackPower = fSignalPower;
    }


    void FakeFreeSpaceSignalGenerator::InitiateEvent(Event* anEvent, int eventID)
    {
        int random_seed_val;
        if ( fRandomSeed != 0 )
        {
            random_seed_val = fRandomSeed;
        }
        else
        {
            std::random_device rd;
            random_seed_val = rd();
        }

        std::geometric_distribution<int> ntracks_distribution(1./fNTracksMean);
        fNTracks = ntracks_distribution(fRandomEngine)+1;

        anEvent->fEventID = eventID;
        anEvent->fLOFrequency = fLO_frequency;
        anEvent->fRandomSeed = random_seed_val;
    }

    bool FakeFreeSpaceSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }



    bool FakeFreeSpaceSignalGenerator::DoGenerateTime( Signal* aSignal )
    {
        InitializeAntennaArray();

    	FileWriter* aRootTreeWriter = RootTreeWriter::get_instance();
    	aRootTreeWriter->SetFilename(fRootFilename);
        aRootTreeWriter->OpenFile("RECREATE");

        const unsigned nChannels = fNChannels;
        const double tLocustStep = 1./aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
        double tTimeOffset = 0;
        const int signalSize = aSignal->TimeSize();
        double tAmplitude, tPhase;
        LMCThreeVector tCurrentElectronPosition;

        for(unsigned eventID = 0; eventID < fNEvents; ++eventID)// event loop.
        {
            Event* anEvent = new Event();
            InitiateEvent(anEvent, eventID);
            Track aTrack;
            SetTrackProperties(aTrack, 0, tTimeOffset);
            anEvent->AddTrack(aTrack);
            
            unsigned tTrackIndex = 0;

            while( tTrackIndex < fNTracks) //loop over tracks in event
            {
                for(unsigned tChannelIndex = 0; tChannelIndex < nChannels; ++tChannelIndex) // over all channels
                {
                	double LO_phase = 0.;
                    double tElectronTime = 0.;
                    unsigned tTrackIndexRange[2] = {static_cast<unsigned>(fStartTime / tLocustStep), static_cast<unsigned>(fEndTime / tLocustStep)};
                    tTrackIndexRange[1] = std::min(tTrackIndexRange[1], aSignal->TimeSize()*aSignal->DecimationFactor());
                    fCurrentFrequency = fStartFrequency;
                    fCurrentRadialPhase = fStartRadialPhase;

                    for( unsigned index = tTrackIndexRange[0]; index < tTrackIndexRange[1]; ++index ) // advance sampling time
                    {

                        //Update times, frequency, electron position, etc.
                        LO_phase += 2.*LMCConst::Pi()*fLO_frequency * tLocustStep;
                        fCurrentFrequency += fSlope * 1.e6/1.e-3 * tLocustStep;
                        fCurrentRadialPhase += 2.*LMCConst::Pi()*fGradBFrequency * tLocustStep;
                        tElectronTime += tLocustStep;
                        tCurrentElectronPosition = GetElectronPosition(fCurrentRadialPhase);

                        tAmplitude = GetCRESAmplitude(tCurrentElectronPosition, tChannelIndex);
                        tPhase = GetCRESPhase(tElectronTime, tCurrentElectronPosition, tChannelIndex);

                        aSignal->LongSignalTimeComplex()[tChannelIndex*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += tAmplitude * cos(tPhase-LO_phase);
                        aSignal->LongSignalTimeComplex()[tChannelIndex*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += tAmplitude * cos(-LMCConst::Pi()/2. + tPhase-LO_phase);

                    }

                } //channel loop

                ++tTrackIndex;
                tTimeOffset = fEndTime;
                SetTrackProperties(aTrack, tTrackIndex, tTimeOffset);

                if(tTrackIndex == fNTracks)
                {
                    break;
                }

                anEvent->AddTrack(aTrack);

            } //track loop

            aRootTreeWriter->WriteEvent(anEvent);
            delete anEvent;
        } //event loop

        aRootTreeWriter->CloseFile();

        return true;
    }

    bool FakeFreeSpaceSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */

/*
 * LMCFakeTrackSignalGenerator.cc
 *
 *  Created on: Aug 8 2018
 *      Author: plslocum
 */

#include <cmath>
#include "LMCFakeTrackSignalGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"
#include <random>
#include <math.h>

using std::string;

namespace locust
{
    LOGGER( lmclog, "FakeTrackSignalGenerator" );

    MT_REGISTER_GENERATOR(FakeTrackSignalGenerator, "fake-track");

    FakeTrackSignalGenerator::FakeTrackSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &FakeTrackSignalGenerator::DoGenerateTime ),
        fSignalPower( 0. ),
        fStartFrequencyMax( 0. ),
        fStartFrequencyMin( 0. ),
        fStartVPhase( 0. ),
        fSlopeMean( 0. ),
        fSlopeStd( 0. ),
        fStartTimeMin( 0. ),
        fStartTimeMax( 0. ),
        fLO_frequency( 0. ),
        fTrackLengthMean( 0. )
    {
        fRequiredSignalState = Signal::kTime;
    }

    FakeTrackSignalGenerator::~FakeTrackSignalGenerator()
    {
    }


    bool FakeTrackSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        if( aParam->has( "signal-power" ) )
            SetSignalPower( aParam->get_value< double >( "signal-power", fSignalPower ) );

        if( aParam->has( "start-frequency-max" ) )
            SetStartFrequencyMax( aParam->get_value< double >( "start-frequency-max", fStartFrequencyMax ) );

        if( aParam->has( "start-frequency-min" ) )
            SetStartFrequencyMin( aParam->get_value< double >( "start-frequency-min", fStartFrequencyMin ) );

        if( aParam->has( "start-vphase" ) )
            SetStartVPhase( aParam->get_value< double >( "start-vphase", fStartVPhase ) );

        if( aParam->has( "slope-mean" ) )
            SetSlopeMean( aParam->get_value< double >( "slope-mean", fSlopeMean ) );

        if( aParam->has( "slope-std" ) )
            SetSlopeStd( aParam->get_value< double >( "slope-std", fSlopeStd ) );

        if( aParam->has( "start-time-max" ) )
            SetStartTimeMax( aParam->get_value< double >( "start-time-max", fStartTimeMax ) );

        if( aParam->has( "start-time-min" ) )
            SetStartTimeMin( aParam->get_value< double >( "start-time-min", fStartTimeMin ) );

        if( aParam->has( "lo-frequency" ) )
            SetFrequency( aParam->get_value< double >( "lo-frequency", fLO_frequency ) );

        if( aParam->has( "track-length-mean" ) )
            SetTrackLengthMean( aParam->get_value< double >( "track-length-mean", fTrackLengthMean ) );



        if( aParam->has( "domain" ) )
        {
            string domain = aParam->get_value( "domain" );
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


    void FakeTrackSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double FakeTrackSignalGenerator::GetSignalPower() const
    {
        return fSignalPower;
    }

    void FakeTrackSignalGenerator::SetSignalPower( double aPower )
    {
        fSignalPower = aPower;
        return;
    }

    double FakeTrackSignalGenerator::GetStartFrequencyMax() const
    {
        return fStartFrequencyMax;
    }

    void FakeTrackSignalGenerator::SetStartFrequencyMax( double aFrequencyMax )
    {
        fStartFrequencyMax = aFrequencyMax;
        return;
    }

    double FakeTrackSignalGenerator::GetStartFrequencyMin() const
    {
        return fStartFrequencyMin;
    }

    void FakeTrackSignalGenerator::SetStartFrequencyMin( double aFrequencyMin )
    {
        fStartFrequencyMin = aFrequencyMin;
        return;
    }

    double FakeTrackSignalGenerator::GetStartVPhase() const
    {
        return fStartVPhase;
    }

    void FakeTrackSignalGenerator::SetStartVPhase( double aPhase )
    {
        fStartVPhase = aPhase;
        return;
    }

    double FakeTrackSignalGenerator::GetSlopeMean() const
    {
        return fSlopeMean;
    }

    void FakeTrackSignalGenerator::SetSlopeMean( double aSlopeMean )
    {
        fSlopeMean = aSlopeMean;
        return;
    }

    double FakeTrackSignalGenerator::GetSlopeStd() const
    {
        return fSlopeStd;
    }

    void FakeTrackSignalGenerator::SetSlopeStd( double aSlopeStd )
    {
        fSlopeStd = aSlopeStd;
        return;
    }

    double FakeTrackSignalGenerator::GetTrackLengthMean() const
    {
        return fTrackLengthMean;
    }

    void FakeTrackSignalGenerator::SetTrackLengthMean( double aTrackLengthMean )
    {
        fTrackLengthMean = aTrackLengthMean;
        return;
    }

    double FakeTrackSignalGenerator::GetStartTimeMin() const
    {
        return fStartTimeMin;
    }

    void FakeTrackSignalGenerator::SetStartTimeMin( double aTimeMin )
    {
        fStartTimeMin = aTimeMin;
        return;
    }

    double FakeTrackSignalGenerator::GetStartTimeMax() const
    {
        return fStartTimeMax;
    }

    void FakeTrackSignalGenerator::SetStartTimeMax( double aTimeMax )
    {
        fStartTimeMax = aTimeMax;
        return;
    }

    double FakeTrackSignalGenerator::GetFrequency() const
    {
        return fLO_frequency;
    }

    void FakeTrackSignalGenerator::SetFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
        return;
    }



    Signal::State FakeTrackSignalGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void FakeTrackSignalGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &FakeTrackSignalGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &FakeTrackSignalGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    bool FakeTrackSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool FakeTrackSignalGenerator::DoGenerateTime( Signal* aSignal )
    {

        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;

        const unsigned nchannels = fNChannels;
        double LO_phase = 0.;
        double dt = 1./aSignal->DecimationFactor()/(RunLengthCalculator1->GetAcquisitionRate()*1.e6);
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::normal_distribution<double> slope_distribution(fSlopeMean,fSlopeStd);
        double slope_val = slope_distribution(generator);
        std::uniform_real_distribution<double> startfreq_distribution(fStartFrequencyMin,fStartFrequencyMax);
        double startfreq_val = startfreq_distribution(generator);
        std::exponential_distribution<double> tracklength_distribution(1./fTrackLengthMean);
        double tracklength_val = tracklength_distribution(generator);
        std::uniform_real_distribution<double> starttime_distribution(fStartTimeMin,fStartTimeMax);
        double starttime_val = starttime_distribution(generator);
        double endtime_val = starttime_val + tracklength_val/sqrt(1+slope_val*slope_val);
        double endfreq_val = startfreq_val + tracklength_val*slope_val/sqrt(1+slope_val*slope_val);
        /* 
        printf("slope_val: %f\n", slope_val);
        printf("tracklength_val: %f\n", tracklength_val);
        printf("starttime_val: %f\n", starttime_val);
        printf("startfreq_val: %f\n", startfreq_val);
        printf("endtime_val: %f\n", endtime_val);
        printf("endfreq_val: %f\n", endfreq_val); getchar();
        */ 
        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
            double voltage_phase = fStartVPhase;
            double track_frequency = startfreq_val;

            for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )

            {
                double time = (double)index/aSignal->DecimationFactor()/(RunLengthCalculator1->GetAcquisitionRate()*1.e6);

                LO_phase += 2.*LMCConst::Pi()*fLO_frequency*dt;

                if ((time > starttime_val) && (time < endtime_val))
                {
                    track_frequency += slope_val*1.e6/1.e-3*dt;
                    voltage_phase += 2.*LMCConst::Pi()*track_frequency*(dt);
 

                    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += sqrt(50.)*sqrt(fSignalPower)*cos(voltage_phase-LO_phase);
                    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += sqrt(50.)*sqrt(fSignalPower)*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);
                }

            }
        }
        delete RunLengthCalculator1;
        return true;
    }

    bool FakeTrackSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */

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


using std::string;

namespace locust
{
    LOGGER( lmclog, "FakeTrackSignalGenerator" );

    MT_REGISTER_GENERATOR(FakeTrackSignalGenerator, "fake-track");

    FakeTrackSignalGenerator::FakeTrackSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &FakeTrackSignalGenerator::DoGenerateTime ),
        fSignalPower( 0. ),
        fStartFrequency( 0. ),
        fStartVPhase( 0. ),
        fSlope( 0. ),
        fStartTime( 0. ),
        fEndTime( 0. ),
        fLO_frequency( 0. )
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

        if( aParam->has( "start-frequency" ) )
            SetStartFrequency( aParam->get_value< double >( "start-frequency", fStartFrequency ) );

        if( aParam->has( "start-vphase" ) )
            SetStartVPhase( aParam->get_value< double >( "start-vphase", fStartVPhase ) );

        if( aParam->has( "slope" ) )
            SetSlope( aParam->get_value< double >( "slope", fSlope ) );

        if( aParam->has( "start-time" ) )
            SetStartTime( aParam->get_value< double >( "start-time", fStartTime ) );

        if( aParam->has( "end-time" ) )
            SetEndTime( aParam->get_value< double >( "end-time", fEndTime ) );

        if( aParam->has( "lo-frequency" ) )
            SetFrequency( aParam->get_value< double >( "lo-frequency", fLO_frequency ) );



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

    double FakeTrackSignalGenerator::GetStartFrequency() const
    {
        return fStartFrequency;
    }

    void FakeTrackSignalGenerator::SetStartFrequency( double aFrequency )
    {
        fStartFrequency = aFrequency;
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

    double FakeTrackSignalGenerator::GetSlope() const
    {
        return fSlope;
    }

    void FakeTrackSignalGenerator::SetSlope( double aSlope )
    {
        fSlope = aSlope;
        return;
    }

    double FakeTrackSignalGenerator::GetStartTime() const
    {
        return fStartTime;
    }

    void FakeTrackSignalGenerator::SetStartTime( double aTime )
    {
        fStartTime = aTime;
        return;
    }

    double FakeTrackSignalGenerator::GetEndTime() const
    {
        return fEndTime;
    }

    void FakeTrackSignalGenerator::SetEndTime( double aTime )
    {
        fEndTime = aTime;
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

        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
            double voltage_phase = fStartVPhase;
            double track_frequency = fStartFrequency;

            for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )

            {
                double time = (double)index/aSignal->DecimationFactor()/(RunLengthCalculator1->GetAcquisitionRate()*1.e6);

                LO_phase += 2.*LMCConst::Pi()*fLO_frequency*dt;

                if ((time > fStartTime) && (time < fEndTime))
                {
                    track_frequency += fSlope*1.e6/1.e-3*dt;
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

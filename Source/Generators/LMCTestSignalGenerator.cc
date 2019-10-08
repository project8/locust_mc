/*
 * LMCTestSignal.cc
 *
 *  Created on: Jan 14 2015
 *      Author: plslocum after nsoblath
 */

#include <cmath>
#include "LMCTestSignalGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "TestSignalGenerator" );

    MT_REGISTER_GENERATOR(TestSignalGenerator, "test-signal");

    TestSignalGenerator::TestSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &TestSignalGenerator::DoGenerateTime ),
        fRF_frequency( 20.1e9 ),
        fAmplitude( 5.e-8 )
    {
        fRequiredSignalState = Signal::kTime;
    }

    TestSignalGenerator::~TestSignalGenerator()
    {
    }


    bool TestSignalGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "rf-frequency" ) )
        {
        SetRFFrequency( aParam.get_value< double >( "rf-frequency", fRF_frequency ) );
        }

        if( aParam.has( "amplitude" ) )
        {
        SetAmplitude( aParam.get_value< double >( "amplitude", fAmplitude ) );
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


    void TestSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double TestSignalGenerator::GetRFFrequency() const
    {
        return fRF_frequency;
    }

    void TestSignalGenerator::SetRFFrequency( double aFrequency )
    {
        fRF_frequency = aFrequency;
        return;
    }

    double TestSignalGenerator::GetAmplitude() const
    {
        return fAmplitude;
    }

    void TestSignalGenerator::SetAmplitude( double aAmplitude )
    {
        fAmplitude = aAmplitude;
        return;
    }



    Signal::State TestSignalGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void TestSignalGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &TestSignalGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &TestSignalGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    bool TestSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool TestSignalGenerator::DoGenerateTime( Signal* aSignal )
    {

        const unsigned nchannels = fNChannels;

        double voltage_phase = 0.;

        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
            for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
            {

                voltage_phase = 2.*LMCConst::Pi()*fRF_frequency*(double)index/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);

                aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += sqrt(50.)*fAmplitude*cos(voltage_phase);
                aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += sqrt(50.)*fAmplitude*cos(-LMCConst::Pi()/2. + voltage_phase);

//printf("signal %d is with acqrate %g, lo %g and rf %g is %g\n", index, fAcquisitionRate, fLO_frequency, fRF_frequency, aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0]); getchar();


            }
        }
        return true;
    }

    bool TestSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */

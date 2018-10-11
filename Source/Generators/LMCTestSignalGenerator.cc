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
        fLO_frequency( 0. ),
        fAmplitude( 0.24 )
    {
        fRequiredSignalState = Signal::kTime;
    }

    TestSignalGenerator::~TestSignalGenerator()
    {
    }


    bool TestSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        SetFrequency( aParam->get_value< double >( "lo-frequency", fLO_frequency ) );
        SetAmplitude( aParam->get_value< double >( "amplitude", fAmplitude ) );


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


    void TestSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double TestSignalGenerator::GetFrequency() const
    {
        return fLO_frequency;
    }

    void TestSignalGenerator::SetFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
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
        fRequiredSignalState = aDomain;  // pls changed == to =.
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

        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;

        const unsigned nchannels = fNChannels;

        double LO_phase = 0.;
        double voltage_phase = 0.;
        double test_frequency = 20.1e9; // Hz

        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
            for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
            {

                LO_phase = 2.*LMCConst::Pi()*fLO_frequency*(double)index/aSignal->DecimationFactor()/(RunLengthCalculator1->GetAcquisitionRate()*1.e6);
                voltage_phase = 2.*LMCConst::Pi()*test_frequency*(double)index/aSignal->DecimationFactor()/(RunLengthCalculator1->GetAcquisitionRate()*1.e6);

                aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += sqrt(50.)*5.e-8*cos(voltage_phase-LO_phase);
                aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += sqrt(50.)*5.e-8*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);


            }
        }
        delete RunLengthCalculator1;
        return true;
    }

    bool TestSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */

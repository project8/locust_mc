/*
 * LMCTestSignal.cc
 *
 *  Created on: Jan 14 2015
 *      Author: plslocum after nsoblath
 */

#include "LMCTestSignalGenerator.hh"

#include "logger.hh"

using std::string;

namespace locust
{
    LOGGER( lmclog, "TestSignalGenerator" );

    MT_REGISTER_GENERATOR(TestSignalGenerator, "test-signal");

    TestSignalGenerator::TestSignalGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &TestSignalGenerator::DoGenerateFreq ),
            fFrequency( 4000. ),
            fAmplitude( 0.24 )
    {
        fRequiredSignalState = Signal::kFreq;
    }

    TestSignalGenerator::~TestSignalGenerator()
    {
    }

    bool TestSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        SetFrequency( aParam->get_value< double >( "frequency", fFrequency ) );
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
        return fFrequency;
    }

    void TestSignalGenerator::SetFrequency( double aFrequency )
    {
        fFrequency = aFrequency;
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


    bool TestSignalGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool TestSignalGenerator::DoGenerateTime( Signal* aSignal ) const
    {
        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;
        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
            aSignal->SignalTime()[index] += fAmplitude*cos(2.*3.1415926*fFrequency*(double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6));  // T=2PI/A
        }
        delete RunLengthCalculator1;
        return true;
    }

    bool TestSignalGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;
        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
        {
            aSignal->SignalFreq()[index][0] += fAmplitude*cos(2.*3.1415926*fFrequency*(double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6));
            aSignal->SignalFreq()[index][1] += fAmplitude*sin(2.*3.1415926*fFrequency*(double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6));
        }
        delete RunLengthCalculator1;
        return true;
    }

} /* namespace locust */

/*
 * LMCTestSignal.cc
 *
 *  Created on: Jan 14 2015
 *      Author: plslocum after nsoblath
 */

#include <cmath>
#include "LMCTestSignalGenerator.hh"
#include "logger.hh"
#include "LMCSimulationController.hh"
#include "LMCConst.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "TestSignalGenerator" );

    MT_REGISTER_GENERATOR(TestSignalGenerator, "test-signal");

    TestSignalGenerator::TestSignalGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &TestSignalGenerator::DoGenerateTime ),
            fFrequency( 4000. ),
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


    bool TestSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool TestSignalGenerator::DoGenerateTime( Signal* aSignal )
    {

        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;


        SimulationController SimulationController1;
        const unsigned nchannels = SimulationController1.GetNChannels();

        double LO_phase = 0.;
        double voltage_phase = 0.;
        double LO_frequency = 20.15e9; // Hz
        double test_frequency = 20.1e9; // Hz

        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
        for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
        {

        	LO_phase = 2.*LMCConst::Pi()*LO_frequency*(double)index/aSignal->DecimationFactor()/(RunLengthCalculator1->GetAcquisitionRate()*1.e6);
            voltage_phase = 2.*LMCConst::Pi()*test_frequency*(double)index/aSignal->DecimationFactor()/(RunLengthCalculator1->GetAcquisitionRate()*1.e6);

        	if (ch==0)
        	{
            aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += 5.e-8*cos(voltage_phase-LO_phase);
            aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += 5.e-8*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);
        	}
        	else
        	{
            aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += 5.e-8*cos(voltage_phase-LO_phase);
            aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += 5.e-8*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);
        	}
//            printf("acq rate is %g\n", RunLengthCalculator1->GetAcquisitionRate()); getchar();
//            printf("array index is %d\n", ch*aSignal->TimeSize() + index);
//            printf("aSignal->SignalTimeComplex()[0] is %g\n", aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][0]); getchar();
        }
        }
        delete RunLengthCalculator1;
        return true;
    }

    bool TestSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;
        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
        {
            aSignal->SignalFreq()[index][0] += fAmplitude*cos(2.*LMCConst::Pi()*fFrequency*(double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6));
            aSignal->SignalFreq()[index][1] += fAmplitude*sin(2.*LMCConst::Pi()*fFrequency*(double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6));
        }
        delete RunLengthCalculator1;
        return true;
    }

} /* namespace locust */

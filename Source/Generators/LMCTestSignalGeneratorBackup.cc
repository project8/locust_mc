/*
 * LMCTestSignalGenerator.cc
 *
 *  Created on: Dec 15, 2014
 *      Author: pslocum
 */

#include "LMCTestSignalGenerator.hh"

#include "LMCLogger.hh"

namespace locust
{

    LMCLOGGER( lmclog, "TestSignalGenerator" );

    MT_REGISTER_GENERATOR(TestSignalGenerator, "test-signal");

    TestSignalGenerator::TestSignalGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &TestSignalGenerator::DoGenerateTime )  //pls:  important
    {
        fRequiredSignalState = Signal::kTime;
    }

    TestSignalGenerator::~TestSignalGenerator()
    {
    }

    bool TestSignalGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;

        return true;
    }

    void TestSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    bool TestSignalGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }



    bool TestSignalGenerator::DoGenerateTime( Signal* aSignal ) const
    {
        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
           aSignal->SignalTime()[index] += 0.24*sin(2.5e4*(double)index/4194304.*0.0209715);  // T=2PI/A
//           if (index<100) printf("TestSignal output at time %g is %f, where TimeSize is %d\n", (double)index/4194304.*0.021, aSignal->SignalTime( index ), aSignal->TimeSize());
        }
        return true;
    }


} /* namespace locust */

/*
 * LMCSimulationController.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#include "LMCSimulationController.hh"

#include "LMCDigitizer.hh"
#include "LMCGenerator.hh"
#include "LMCLogger.hh"
#include "LMCParam.hh"

namespace locust
{
    LMCLOGGER( lmclog, "SimulationController" );

    SimulationController::SimulationController() :
            fRNG(),
            fFirstGenerator( NULL ),
            fRunLengthCalc(),
            fEggWriter()
    {
        SetRNGSeed();
    }

    SimulationController::~SimulationController()
    {
    }

    bool SimulationController::Configure( const ParamNode* aNode )
    {
        if( aNode == NULL ) return true;

        if( aNode->Has( "rng-seed" ) )
        {
            SetRNGSeed( aNode->GetValue< int >( "rng-seed" ) );
        }

        // configure the run-length calculator
        if( ! fRunLengthCalc.Configure( aNode ) )
        {
            LMCERROR( lmclog, "Error configuring the run length calculator" );
            return false;
        }

        // configure the egg writer
        if( ! fEggWriter.Configure( aNode ) )
        {
            LMCERROR( lmclog, "Error configuring the egg writer" );
            return false;
        }

        return true;
    }

    void SimulationController::SetRNGSeed(int seed)
    {
        fRNG.Reseed( seed );
        return;
    }

    void SimulationController::SetRNGSeed()
    {
        fRNG.Reseed( RandomLib::Random::SeedWord() );
        return;
    }


    void SimulationController::SetFirstGenerator( Generator* firstGen )
    {
        fFirstGenerator = firstGen;
        fRunLengthCalc.SetFirstGenerator( fFirstGenerator );
        return;
    }

    bool SimulationController::Prepare()
    {
        LMCINFO( lmclog, "Preparing for run" );

        Generator* nextGenerator = fFirstGenerator;
        while( nextGenerator != NULL )
        {
            nextGenerator->SetRNG( &fRNG );
            nextGenerator = nextGenerator->GetNextGenerator();
        }

        if( ! fRunLengthCalc.VisitGenerators() )
        {
            LMCERROR( lmclog, "Error while the run-length calculator was visiting generators" );
            return false;
        }

        // do the final determination of the run length
        if( ! fRunLengthCalc.CalculateRunLength() )
        {
            LMCERROR( lmclog, "Error while the run-length calculator was calculating the run length" );
            return false;
        }

        // prepare the egg file (writes header, allocates record memory, etc)
        if( ! fEggWriter.PrepareEgg( &fRunLengthCalc, FindDigitizer() ) )
        {
            LMCERROR( lmclog, "Error preparing the egg file" );
            return false;
        }

        return true;
    }

    bool SimulationController::Run()
    {
        if( fFirstGenerator == NULL )
        {
            LMCWARN( lmclog, "First generator is not present" );
            return false;
        }

        unsigned nRecords = fRunLengthCalc.GetNRecords();
        unsigned recordSize = fRunLengthCalc.GetRecordSize();

        LMCINFO( lmclog, "Commencing the run" );

        for( unsigned record = 0; record < nRecords; ++record )
        {
            LMCINFO( lmclog, "Simulating record " << record );
            Signal* simulatedSignal = fFirstGenerator->Run( recordSize );
            if( simulatedSignal == NULL )
            {
                LMCERROR( lmclog, "Signal was not simulated" );
                return false;
            }


            if( simulatedSignal->GetDigitalIsSigned() )
            {
                for( unsigned index = 0; index < 10; ++index )  // pls changed loop range to 10.
                {
                    LMCWARN( lmclog, index << "  " << simulatedSignal->SignalTime()[index] << "  " << (int)simulatedSignal->SignalDigitalS()[index] );  // pls added (int)
                }
            }
            else
            {
                for( unsigned index = 0; index < 10; ++index )  // pls changed loop range to 10.
                {
                    LMCWARN( lmclog, index << "  " << simulatedSignal->SignalTime()[index] << "  " << (int)simulatedSignal->SignalDigitalUS()[index] );  // pls added (int)
                }
            }

            if( ! fEggWriter.WriteRecord( simulatedSignal ) )
            {
                LMCERROR( lmclog, "Something went wrong while writing record " << record );
                delete simulatedSignal;
                return false;
            }

            // temporarily, immediately cleanup
            delete simulatedSignal;
        }

        return true;
    }

    bool SimulationController::Finalize()
    {
        LMCINFO( lmclog, "Finalizing the run" );

        if(! fEggWriter.FinalizeEgg() )
        {
            LMCERROR( lmclog, "Error while finalizing the egg file" );
            return false;
        }

        return true;
    }

    const Digitizer* SimulationController::FindDigitizer() const
    {
        const Generator* thisGenerator = fFirstGenerator;
        const Digitizer* asDigitizer = dynamic_cast< const Digitizer* >( thisGenerator );
        while( thisGenerator != NULL && asDigitizer == NULL )
        {
            thisGenerator = thisGenerator->GetNextGenerator();
            asDigitizer = dynamic_cast< const Digitizer* >( thisGenerator );
        }
        return asDigitizer;
    }

} /* namespace locust */

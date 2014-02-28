/*
 * LMCSimulationController.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#include "LMCSimulationController.hh"

#include "LMCGenerator.hh"
#include "LMCLogger.hh"
#include "LMCParam.hh"

namespace locust
{
    LMCLOGGER( lmclog, "SimulationController" );

    SimulationController::SimulationController() :
            fFirstGenerator( NULL ),
            fRunLengthCalc(),
            fEggWriter()
    {
    }

    SimulationController::~SimulationController()
    {
    }

    bool SimulationController::Configure( const ParamNode* aNode )
    {
        if( aNode == NULL ) return true;

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

    void SimulationController::SetFirstGenerator( const Generator* firstGen )
    {
        fFirstGenerator = firstGen;
        return;
    }

    bool SimulationController::Prepare()
    {
        LMCINFO( lmclog, "Preparing for run" );

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

        // overwrite the run duration in the egg writer, then prepare the egg file
        fEggWriter.SetDuration( fRunLengthCalc.GetDuration() );
        if( ! fEggWriter.PrepareEgg() )
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
            for( unsigned index = 0; index < 100; ++index )
            {
                LMCWARN( lmclog, simulatedSignal->SignalTime( index ) );
            }

            if( ! simulatedSignal->ToState( Signal::kDigital ) )
            {
                LMCERROR( lmclog, "Please digitize the signal before attempting to write to an egg file" );
                delete simulatedSignal;
                return false;
            }
            if( ! fEggWriter.WriteRecord( simulatedSignal ) )
            {
                LMCERROR( lmclog, "Something went wrong while writing record " << index );
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

} /* namespace locust */

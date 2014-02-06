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
            fNRecords( 1 ),
            fDuration( 0.021 ),
            fRecordSize( 4194304 ),
            fBinWidth( 5.e9 ),
            fFirstGenerator( NULL )
    {
    }

    SimulationController::~SimulationController()
    {
    }

    void SimulationController::Configure( const ParamNode* aNode )
    {
        if( aNode == NULL ) return;

        // NOTE: fBinWidth, fRecordSize, and fNRecords are set before setting the duration on purpose.
        // The former three all modify the duration.

        if( aNode->has( "digitizer-rate" ) )
            SetBinWidth( 1. / aNode->get_value< double >( "digitizer-rate" ) );
        else if( aNode->has( "bin-width" ) )
            SetBinWidth( aNode->get_value< double >( "bin-width" ) );

        if( aNode->has( "record-size" ) )
            SetRecordSize( aNode->get_value< unsigned >( "record-size" ) );

        if( aNode->has( "n-records" ) )
            SetNRecords( aNode->get_value< unsigned >( "n-records" ) );

        if( aNode->has( "duration" ) )
            SetDuration( aNode->get_value< double >( "duration" ) );

        return;
    }

    unsigned SimulationController::GetNRecords() const
    {
        return fNRecords;
    }

    void SimulationController::SetNRecords( unsigned recs )
    {
        fNRecords = recs;
        fDuration = fNRecords * fBinWidth * fRecordSize;
        return;
    }

    double SimulationController::GetDuration() const
    {
        return fDuration;
    }

    void SimulationController::SetDuration( double duration )
    {
        fNRecords = duration / double(fBinWidth * fRecordSize);
        fDuration = fNRecords * fBinWidth * fRecordSize;
        return;
    }

    unsigned SimulationController::GetRecordSize() const
    {
        return fRecordSize;
    }

    void SimulationController::SetRecordSize( unsigned size )
    {
        fRecordSize = size;
        fDuration = fNRecords * fBinWidth * fRecordSize;
        return;
    }

    double SimulationController::GetBinWidth() const
    {
        return fBinWidth;
    }

    void SimulationController::SetBinWidth( double bw )
    {
        fBinWidth = bw;
        fDuration = fNRecords * fBinWidth * fRecordSize;
        return;
    }

    void SimulationController::SetFirstGenerator( const Generator* firstGen )
    {
        fFirstGenerator = firstGen;
        return;
    }

    void SimulationController::Run() const
    {
        if( fFirstGenerator == NULL )
        {
            LMCWARN( lmclog, "First generator is not present" );
            return;
        }

        for( unsigned record = 0; record < fNRecords; ++record )
        {
            LMCINFO( lmclog, "Simulating record " << record );
            Signal* simulatedSignal = fFirstGenerator->Run( fRecordSize );
            for( unsigned index = 0; index < 100; ++index )
            {
                LMCWARN( lmclog, simulatedSignal->SignalTime( index ) );
            }
            // temporarily, immediately cleanup
            delete simulatedSignal;
        }

        return;
    }



} /* namespace locust */

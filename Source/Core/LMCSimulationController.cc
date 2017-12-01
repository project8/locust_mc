/*
 * LMCSimulationController.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#include "LMCSimulationController.hh"

#include "LMCDigitizer.hh"
#include "LMCGenerator.hh"
#include "logger.hh"
#include "param.hh"

namespace locust
{
    LOGGER( lmclog, "SimulationController" );

    SimulationController::SimulationController() :
            fFirstGenerator( NULL ),
            fRunLengthCalc(),
            fEggWriter()
    {
        SetRNGSeed();
    }

    SimulationController::~SimulationController()
    {
    }

    bool SimulationController::Configure( const scarab::param_node* aNode )
    {
        if( aNode == NULL ) return true;

        if( aNode->has( "rng-seed" ) )
        {
            SetRNGSeed( aNode->get_value< int >( "rng-seed" ) );
        }

        // configure the run-length calculator
        if( ! fRunLengthCalc.Configure( aNode ) )
        {
            LERROR( lmclog, "Error configuring the run length calculator" );
            return false;
        }

        // configure the egg writer
        if( ! fEggWriter.Configure( aNode ) )
        {
            LERROR( lmclog, "Error configuring the egg writer" );
            return false;
        }

        return true;
    }

    void SimulationController::SetRNGSeed(int seed)
    {
        Generator::RNG().seed( seed );
        return;
    }

    void SimulationController::SetRNGSeed()
    {
        Generator::RNG().seed();
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
        LINFO( lmclog, "Preparing for run" );

        if( ! fRunLengthCalc.VisitGenerators() )
        {
            LERROR( lmclog, "Error while the run-length calculator was visiting generators" );
            return false;
        }

        // do the final determination of the run length
        if( ! fRunLengthCalc.CalculateRunLength() )
        {
            LERROR( lmclog, "Error while the run-length calculator was calculating the run length" );
            return false;
        }

        // prepare the egg file (writes header, allocates record memory, etc)
        if( ! fEggWriter.PrepareEgg( &fRunLengthCalc, FindDigitizer() ) )
        {
            LERROR( lmclog, "Error preparing the egg file" );
            return false;
        }

        return true;
    }

    bool SimulationController::Run()
    {
        if( fFirstGenerator == NULL )
        {
            LWARN( lmclog, "First generator is not present" );
            return false;
        }

        bool IQStream = true;  // make this a parameter in the json file.
        unsigned nRecords = fRunLengthCalc.GetNRecords();
        unsigned recordSize = fRunLengthCalc.GetRecordSize();

        LINFO( lmclog, "Commencing the run" );

        for( unsigned record = 0; record < nRecords; ++record )
        {
            LINFO( lmclog, "Simulating record " << record );
            Signal* simulatedSignal = fFirstGenerator->Run( recordSize );
            if( simulatedSignal == NULL )
            {
                LERROR( lmclog, "Signal was not simulated" );
                return false;
            }



            if( simulatedSignal->GetDigitalIsSigned() )
            {
                for( unsigned index = 0; index < 10; ++index )  // pls changed loop range to 10.
                {
                	if (!IQStream)
                	{
                        LWARN( lmclog, index << "  " << simulatedSignal->SignalTime()[index] << "  " << (int)simulatedSignal->SignalDigitalS()[index] );  // pls added (int)
                	}
                	else
                	{
                        LWARN( lmclog, index << "  " << simulatedSignal->SignalTimeComplex()[index][0] << "  " << (int)simulatedSignal->SignalDigitalS()[index*2] );  // pls added (int)
                        LWARN( lmclog, index << "  " << simulatedSignal->SignalTimeComplex()[index][1] << "  " << (int)simulatedSignal->SignalDigitalS()[index*2+1] );  // pls added (int)
                	}
                }
            }
            else
            {
                for( unsigned index = 0; index < 10; ++index )  // pls changed loop range to 10.
                {
                	if (!IQStream)
                	{
                        LWARN( lmclog, index << "  " << simulatedSignal->SignalTime()[index] << "  " << (int)simulatedSignal->SignalDigitalUS()[index] );  // pls added (int)
                	}
                	else
                	{
                        LWARN( lmclog, index << "  " << simulatedSignal->SignalTimeComplex()[index][0] << "  " << (int)simulatedSignal->SignalDigitalUS()[index*2] );  // pls added (int)
                        LWARN( lmclog, index << "  " << simulatedSignal->SignalTimeComplex()[index][1] << "  " << (int)simulatedSignal->SignalDigitalUS()[index*2+1] );  // pls added (int)
                	}
                }
            }

            if( ! fEggWriter.WriteRecord( simulatedSignal ) )
            {
                LERROR( lmclog, "Something went wrong while writing record " << record );
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
        LINFO( lmclog, "Finalizing the run" );

        if(! fEggWriter.FinalizeEgg() )
        {
            LERROR( lmclog, "Error while finalizing the egg file" );
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

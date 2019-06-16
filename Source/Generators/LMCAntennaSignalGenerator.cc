/*
 * LMCAntennaSignalGenerator.cc
 *
 *  Created on: May 4, 2019
 *      Author: pslocum
 */

#include "LMCAntennaSignalGenerator.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <math.h>       /* sin */
#include "LMCGlobalsDeclaration.hh"
#include "LMCDigitizer.hh"

namespace locust
{
    LOGGER( lmclog, "AntennaSignalGenerator" );

    MT_REGISTER_GENERATOR(AntennaSignalGenerator, "antenna-signal");

    AntennaSignalGenerator::AntennaSignalGenerator( const std::string& aName ) :
        Generator( aName ),
	fInputSignalType(1),
	fInputFrequency(27.0), //Assume 27 
        fInputAmplitude(1)
    {
        fRequiredSignalState = Signal::kTime;
    }

    AntennaSignalGenerator::~AntennaSignalGenerator()
    {
    }

    bool AntennaSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;
	 
	if(!fFieldEstimator.Configure(aParam))
	{
		LERROR(lmclog,"Error configuring field estimator class");
	}

	if( aParam->has( "input-signal-type" ) )
        {
            fInputSignalType = aParam->get_value< int >( "input-signal-type" );
        }

	if( aParam->has( "input-signal-frequency" ) )
        {
            fInputSignalType = aParam->get_value< double >( "input-signal-frequency" );
        }

	if( aParam->has( "input-signal-amplitude" ) )
        {
            fInputSignalType = aParam->get_value< double >( "input-signal-amplitude" );
        }

        return true;
    }

    void AntennaSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    void AntennaSignalGenerator::GenerateSignal(Signal* aSignal)
    {
	if(fInputSignalType==1) // sin wave
	{
	     for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
	     {
		 aSignal->LongSignalTimeComplex()[index][0] = fInputAmplitude* sin(index);
		 aSignal->LongSignalTimeComplex()[index][1] = fInputAmplitude* cos(index);
	     }
	}


	else //Else case also sin for now still
	{
	     for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
	     {
		 aSignal->LongSignalTimeComplex()[index][0] = fInputAmplitude* sin(index);
		 aSignal->LongSignalTimeComplex()[index][1] = fInputAmplitude* cos(index);
	     }
	}
    }

    bool AntennaSignalGenerator::DoGenerate( Signal* aSignal )
    {
	fFieldEstimator.ReadFIRFile();
	GenerateSignal(aSignal);

        return true;
    }

} /* namespace locust */

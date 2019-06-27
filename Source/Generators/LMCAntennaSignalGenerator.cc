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

using std::string;

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
            fInputFrequency= aParam->get_value< double >( "input-signal-frequency" );
        }

	if( aParam->has( "input-signal-amplitude" ) )
        {
            fInputAmplitude = aParam->get_value< double >( "input-signal-amplitude" );
        }

	// PTS: Need to check if this needs to be set separately for both generators as well as the transmitter
	// Implement SetDomain functionality later
	//Defining the domain to be used
	/*
        if( aParam->has( "domain" ) )
	{
	    string domain = aParam->get_value( "domain" );
            if( domain == "time" )
	    {
		SetDomain( Signal::kTime );
		LDEBUG( lmclog, "Generating simulation in time domain.");
	    }
            else if( domain == "freq" )
	    {
		SetDomain( Signal::kFreq );
		LDEBUG( lmclog, "Generating simulation in frequency domain.");
	    }
            else
	    {
		LDEBUG( lmclog, "Unable to use domain requested: <" << domain << ">");
		return false;
	    }

	}*/
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

	else //Else case also sin for now
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
	InitializeBuffers(fFieldEstimator.GetFilterSize());
	GenerateSignal(aSignal);	
	fFieldEstimator.ConvolveWithFIRFilter(aSignal);
        return true;
    }

    void AntennaSignalGenerator::InitializeBuffers(unsigned filterbuffersize)
    {
    	FieldBuffer aFieldBuffer;
	delayedVoltageBuffer = aFieldBuffer.InitializeBuffer(1,1,filterbuffersize);
	//const unsigned nchannels = fNChannels;
    }
} /* namespace locust */

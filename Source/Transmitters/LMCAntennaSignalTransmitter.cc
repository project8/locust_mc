/*
 * LMCAntennaSignalTransmitter.cc
 *
 *  Created on: May 4, 2019
 *      Author: pslocum
 */

#include "LMCAntennaSignalTransmitter.hh"

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
    LOGGER( lmclog, "AntennaSignalTransmitter" );

    AntennaSignalTransmitter::AntennaSignalTransmitter() :
	fInputSignalType(1),
	fInputFrequency(27.0e9), //Should use this for the rf frequency ?
        fInputAmplitude(1)
    {
    }

    AntennaSignalTransmitter::~AntennaSignalTransmitter()
    {
    }

    bool AntennaSignalTransmitter::Configure( const scarab::param_node* aParam )
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
	return true;
    }

    double AntennaSignalTransmitter::GenerateSignal(Signal *aSignal,double acquisitionRate)
    {
	double estimatedField=0.0;
	double voltagePhase=fPhaseDelay;
	if(fInputSignalType==1) // sin wave
	{
	     for( unsigned index = 0; index <fFieldEstimator.GetFilterSize();index++)
	     {
		 delayedVoltageBuffer[0].push_back(fInputAmplitude* cos(voltagePhase));
	     	 delayedVoltageBuffer[0].pop_front();
		 voltagePhase += 2.*LMCConst::Pi()*fInputFrequency*fFieldEstimator.GetFilterdt();
	     }
	}

	else //Else case also sin for now
	{
	}
	estimatedField=fFieldEstimator.ConvolveWithFIRFilter(delayedVoltageBuffer[0]);
	//std::cout<< "1: "<<fPhaseDelay <<std::endl;
	fPhaseDelay+= 2.*LMCConst::Pi()*fInputFrequency/aSignal->DecimationFactor()/(acquisitionRate*1.e6);
	return estimatedField;
    }

    bool AntennaSignalTransmitter::InitializeTransmitter()
    {
	fFieldEstimator.ReadFIRFile();
	double filterSize=fFieldEstimator.GetFilterSize();
	InitializeBuffers(filterSize);
	//This has to be added later
	fInitialPhaseDelay = -2.*LMCConst::Pi()*filterSize*fFieldEstimator.GetFilterdt()*fInputFrequency;
	fPhaseDelay = fInitialPhaseDelay;
	return true;
    }

    double AntennaSignalTransmitter::GetInitialPhaseDelay()
    {
	    return fInitialPhaseDelay;
    }

    void AntennaSignalTransmitter::InitializeBuffers(unsigned filterbuffersize)
    {
    	FieldBuffer aFieldBuffer;
	delayedVoltageBuffer = aFieldBuffer.InitializeBuffer(1,1,filterbuffersize);
	//const unsigned nchannels = fNChannels;
    }
} /* namespace locust */

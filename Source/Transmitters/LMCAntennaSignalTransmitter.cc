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
	fInputFrequency(27.0e9), //Assume 27 
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

    double AntennaSignalTransmitter::GenerateSignal(Signal *aSignal,double acquisitionRate)
    {
	double estimatedField=0.0;
	double voltagePhase=0.0;
	//std::cout<< "1: "<<phaseDelay<<" : "<<voltagePhase<<" : " << timeNumber<< std::endl;
	if(fInputSignalType==1) // sin wave
	{
	     for( unsigned index = 0; index <fFieldEstimator.GetFilterSize();index++)
	     {
	     	 voltagePhase += 2.*LMCConst::Pi()*fInputFrequency*fFieldEstimator.GetFilterdt()+phaseDelay;
		 delayedVoltageBuffer[0].push_back(fInputAmplitude* cos(voltagePhase));
	     	 delayedVoltageBuffer[0].pop_front();
	     }
	}

	else //Else case also sin for now
	{
	}
	//std::cout<< "2: "<<phaseDelay<<" : "<<voltagePhase<<" : " << timeNumber<< std::endl;
	estimatedField=fFieldEstimator.ConvolveWithFIRFilter(aSignal);
	phaseDelay+= 2.*LMCConst::Pi()*fInputFrequency/aSignal->DecimationFactor()/(acquisitionRate*1.e6);
	timeNumber+=1;
	//std::cout<< "3: "<<phaseDelay<<" : "<<voltagePhase<<" : " << timeNumber<< std::endl;
	return estimatedField;
    }

    bool AntennaSignalTransmitter::InitializeTransmitter()
    {
	fFieldEstimator.ReadFIRFile();
	double filterSize=fFieldEstimator.GetFilterSize();
	InitializeBuffers(filterSize);
	phaseDelay= -2.*LMCConst::Pi()*filterSize*fFieldEstimator.GetFilterdt()*fInputFrequency;
	std::cout <<filterSize<<": "<<fFieldEstimator.GetFilterdt()<<":"<<fInputFrequency <<"  "<<phaseDelay/2/LMCConst::Pi()<<std::endl;
	return true;
    }

    void AntennaSignalTransmitter::InitializeBuffers(unsigned filterbuffersize)
    {
    	FieldBuffer aFieldBuffer;
	delayedVoltageBuffer = aFieldBuffer.InitializeBuffer(1,1,filterbuffersize);
	//const unsigned nchannels = fNChannels;
    }
} /* namespace locust */

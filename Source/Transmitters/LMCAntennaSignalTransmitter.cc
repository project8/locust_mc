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
#include <math.h>       
#include "LMCGlobalsDeclaration.hh"
#include "LMCDigitizer.hh"

using std::string;

namespace locust
{
    LOGGER( lmclog, "AntennaSignalTransmitter" );

    AntennaSignalTransmitter::AntennaSignalTransmitter() :
	fInputSignalType(1),
	fInputFrequency(27.0e9), //Should be the samne as the value used in the dipole signal generator
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

	if( aParam->has( "array-radius" ) )
	{
		fArrayRadius = aParam->get_value< double >( "array-radius" );
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
	if(fInputSignalType==1) //sinusoidal wave for dipole antenna
	{
	     for( unsigned index = 0; index <fFieldEstimator.GetFilterSize();index++)
	     {
		 double voltageValue = fFieldEstimator.GetFieldAtOrigin(fInputAmplitude,voltagePhase);
		 delayedVoltageBuffer[0].push_back(voltageValue);
	     	 delayedVoltageBuffer[0].pop_front();
		 voltagePhase += 2.*LMCConst::Pi()*fInputFrequency*fFieldEstimator.GetFilterResolution();
	     }
	}

	else// For now using sinusoidal as well 
	{
	     for( unsigned index = 0; index <fFieldEstimator.GetFilterSize();index++)
	     {
		 double voltageValue = fFieldEstimator.GetFieldAtOrigin(fInputAmplitude,voltagePhase);
		 delayedVoltageBuffer[0].push_back(voltageValue);
	     	 delayedVoltageBuffer[0].pop_front();
		 voltagePhase += 2.*LMCConst::Pi()*fInputFrequency*fFieldEstimator.GetFilterResolution();
	     }
	}
	estimatedField=fFieldEstimator.ConvolveWithFIRFilter(delayedVoltageBuffer[0]);
	fPhaseDelay+= 2.*LMCConst::Pi()*fInputFrequency/aSignal->DecimationFactor()/(acquisitionRate*1.e6);
	return estimatedField;
    }

    bool AntennaSignalTransmitter::InitializeTransmitter()
    {
	if(!fFieldEstimator.ReadFIRFile())
	{
		return false;
	}
	double filterSize=fFieldEstimator.GetFilterSize();
	InitializeBuffers(filterSize);
	fInitialPhaseDelay = -2.*LMCConst::Pi()*(filterSize*fFieldEstimator.GetFilterResolution()+fArrayRadius/LMCConst::C())*fInputFrequency;
	fPhaseDelay = fInitialPhaseDelay+travelPhaseDelay;
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
    }
} /* namespace locust */

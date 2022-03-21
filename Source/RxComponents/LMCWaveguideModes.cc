/*
 * LMCWaveguideModes.cc
 *
 *  Created on: Mar 17, 2022
 *      Author: pslocum
 */

#include "LMCWaveguideModes.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "WaveguideModes" );

    WaveguideModes::WaveguideModes()
    {
    }

    WaveguideModes::~WaveguideModes()
    {
    }


    bool WaveguideModes::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombiner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from WaveguideModes subclass");
    		return false;
    	}

    	return true;
    }

	bool WaveguideModes::AddOneModeToCavityProbe(Signal* aSignal, double excitationAmplitude, double EFieldAtProbe, double dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex)
	{

		fVoltagePhase += dopplerFrequency * dt;
		double voltageValue = excitationAmplitude * cos(fVoltagePhase);

		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);


		if ( GetVoltageCheck() && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );

		return true;
	}


} /* namespace locust */

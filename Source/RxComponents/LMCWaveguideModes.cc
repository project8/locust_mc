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

    WaveguideModes::WaveguideModes():
		fVoltagePhaseAntenna( 0. ),
		fVoltagePhaseShort( 0. )
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

    double WaveguideModes::GroupVelocity(double fcyc, double aDimX)
    {
    	// fcyc in cycles/sec
        double cutOffFrequency = LMCConst::C() * LMCConst::Pi() / aDimX; // a in m
        double groupVelocity = LMCConst::C() * sqrt( 1. - pow(cutOffFrequency/( 2.*LMCConst::Pi()*fcyc  ), 2.) );
    	return groupVelocity;
    }

    bool WaveguideModes::InitializeVoltagePhases(std::vector<double> tKassParticleXP, std::vector<double> dopplerFrequency, double aCenterToAntenna, double aCenterToShort, double aDimX)
    {
    	double fcyc = tKassParticleXP[7]/2./LMCConst::Pi(); // cycles/sec
    	double tPositionZ = tKassParticleXP[2];

        double tVoltagePhaseAntenna = 2.*LMCConst::Pi()*(aCenterToAntenna - tPositionZ) / (GroupVelocity(fcyc, aDimX) / dopplerFrequency[0]);
        SetVoltagePhaseAntenna( tVoltagePhaseAntenna );

        double tVoltagePhaseShort = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(aCenterToShort + aCenterToAntenna) /
                (GroupVelocity(fcyc, aDimX) / dopplerFrequency[1]);  // phase of reflected field at antenna.
        SetVoltagePhaseShort( tVoltagePhaseShort );

    	return true;
    }

	bool WaveguideModes::AddOneModeToCavityProbe(Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex)
	{

		double dopplerFrequencyAntenna = dopplerFrequency[0];
		double dopplerFrequencyShort = dopplerFrequency[1];

		SetVoltagePhaseAntenna( GetVoltagePhaseAntenna() + dopplerFrequencyAntenna * dt);
		SetVoltagePhaseShort( GetVoltagePhaseShort() + dopplerFrequencyShort * dt);

		double voltageValue = excitationAmplitude;

		if ( GetWaveguideShortIsPresent() ) // with short:
		{
			voltageValue *= ( cos(GetVoltagePhaseAntenna()) + cos(GetVoltagePhaseShort()) );
			aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
	    	aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);
		}
		else // without short:
		{
			voltageValue *= cos(GetVoltagePhaseAntenna());
    		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
	    	aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);
		}

		if ( GetVoltageCheck() && (sampleIndex%1000 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );

		return true;
	}

    double WaveguideModes::GetVoltagePhaseAntenna()
    {
    	return fVoltagePhaseAntenna;
    }
    void WaveguideModes::SetVoltagePhaseAntenna ( double aPhase )
    {
        fVoltagePhaseAntenna = aPhase;
    }

    double WaveguideModes::GetVoltagePhaseShort()
    {
    	return fVoltagePhaseShort;
    }
    void WaveguideModes::SetVoltagePhaseShort ( double aPhase )
    {
        fVoltagePhaseShort = aPhase;
    }



} /* namespace locust */

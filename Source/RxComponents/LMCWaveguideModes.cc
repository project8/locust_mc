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

    double WaveguideModes::GroupVelocity(double fcyc, double aDimX)
    {
    	// fcyc in cycles/sec
        double cutOffFrequency = LMCConst::C() * LMCConst::Pi() / aDimX; // a in m
        double groupVelocity = LMCConst::C() * sqrt( 1. - pow(cutOffFrequency/( 2.*LMCConst::Pi()*fcyc  ), 2.) );
    	return groupVelocity;
    }

    bool WaveguideModes::InitializeVoltagePhases(std::vector<double> tKassParticleXP, double aDopplerFrequencyAntenna, double aDopplerFrequencyShort, double aCenterToAntenna, double aCenterToShort, double aDimX)
    {
    	double fcyc = tKassParticleXP[7]/2./LMCConst::Pi(); // cycles/sec
    	double tPositionZ = tKassParticleXP[2];

        fVoltagePhaseAntenna = 2.*LMCConst::Pi()*(aCenterToAntenna - tPositionZ) / (GroupVelocity(fcyc, aDimX) / aDopplerFrequencyAntenna);

        fVoltagePhaseShort = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(aCenterToShort + aCenterToAntenna) /
                (GroupVelocity(fcyc, aDimX) / aDopplerFrequencyShort);  // phase of reflected field at antenna.

    	return true;
    }

	bool WaveguideModes::AddOneModeToCavityProbe(Signal* aSignal, double excitationAmplitude, double EFieldAtProbe, double dopplerFrequencyAntenna, double dopplerFrequencyShort, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, double eventTime)
	{

		fVoltagePhaseAntenna += dopplerFrequencyAntenna * dt;
		fVoltagePhaseShort += dopplerFrequencyShort * dt;

		double voltageValue = excitationAmplitude * (cos(fVoltagePhaseAntenna) + cos(fVoltagePhaseShort));

		if ( !GetWaveguideShortIsPresent() ) // no short:
		{
			// override default case, omitting reflected signal:
			voltageValue = excitationAmplitude * cos(fVoltagePhaseAntenna);
		}

		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);


		if ( GetVoltageCheck() && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );

		return true;
	}


} /* namespace locust */

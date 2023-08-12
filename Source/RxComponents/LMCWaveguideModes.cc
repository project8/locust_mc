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
		fVoltagePhaseShort( 0. ),
	    fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
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

    	SizeNChannels(GetNChannels());

    	return true;
    }

    bool WaveguideModes::SizeNChannels(int aNumberOfChannels)
    {

    	SetNChannels(aNumberOfChannels);

    	std::vector<std::vector<std::vector<std::vector<double>>>> tZeroVector;
    	fVoltagePhaseAntenna.swap(tZeroVector);
    	fVoltagePhaseShort.swap(tZeroVector);
    	fVoltagePhaseAntenna.resize(aNumberOfChannels);
    	fVoltagePhaseShort.resize(aNumberOfChannels);

    	for (int n = 0; n < GetNChannels(); n++)
    	{
    		fVoltagePhaseAntenna[n].resize(GetNCavityModes());
    		fVoltagePhaseShort[n].resize(GetNCavityModes());
    		for (int i = 0; i < GetNCavityModes(); i++)
    		{
    			fVoltagePhaseAntenna[n][i].resize(GetNCavityModes());
    			fVoltagePhaseShort[n][i].resize(GetNCavityModes());
    			for (int j = 0; j < GetNCavityModes(); j++)
    			{
    				fVoltagePhaseAntenna[n][i][j].resize(GetNCavityModes());
    				fVoltagePhaseShort[n][i][j].resize(GetNCavityModes());
    			}
    		}
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

    bool WaveguideModes::InitializeVoltagePhases(std::vector<double> tKassParticleXP, std::vector<double> dopplerFrequency, int aChannel, int l, int m, int n)
    {
    	double fcyc = tKassParticleXP[7]/2./LMCConst::Pi(); // cycles/sec
    	double tPositionZ = tKassParticleXP[2];
    	double aCenterToAntenna = fInterface->fCENTER_TO_ANTENNA;
    	double aCenterToShort = fInterface->fCENTER_TO_SHORT;
    	double aDimX = fInterface->fField->GetDimX();

        double tVoltagePhaseAntenna = 2.*LMCConst::Pi()*(aCenterToAntenna - tPositionZ) / (GroupVelocity(fcyc, aDimX) / dopplerFrequency[0]);
        SetVoltagePhaseAntenna( tVoltagePhaseAntenna, aChannel, l, m, n );

        double tVoltagePhaseShort = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(aCenterToShort + aCenterToAntenna) /
                (GroupVelocity(fcyc, aDimX) / dopplerFrequency[1]);  // phase of reflected field at antenna.
        SetVoltagePhaseShort( tVoltagePhaseShort, aChannel, l, m, n );

    	return true;
    }

	bool WaveguideModes::AddOneModeToCavityProbe(int l, int m, int n, Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, int channelIndex, bool initParticle)
	{

		double dopplerFrequencyAntenna = dopplerFrequency.front();
		double dopplerFrequencyShort = dopplerFrequency.back();

		if (initParticle)
		{
			InitializeVoltagePhases(particleXP, dopplerFrequency, channelIndex, l, m, n);
		}

		double newPhaseAntenna = GetVoltagePhaseAntenna(channelIndex, l, m, n) + dopplerFrequencyAntenna * dt;
		double newPhaseShort = GetVoltagePhaseShort(channelIndex, l, m, n) + dopplerFrequencyShort * dt;
		SetVoltagePhaseAntenna( newPhaseAntenna, channelIndex, l, m, n);
		SetVoltagePhaseShort( newPhaseShort, channelIndex, l, m, n);

		double voltageValue = excitationAmplitude;

		if ( GetWaveguideShortIsPresent() ) // with short:
		{
			voltageValue *= ( cos(GetVoltagePhaseAntenna(channelIndex, l, m, n)) + cos(GetVoltagePhaseShort(channelIndex, l, m, n)) );
			aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
	    	aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);
		}
		else // without short:
		{
			voltageValue *= cos(GetVoltagePhaseAntenna(channelIndex, l, m, n));
    		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
	    	aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);
		}

		if ( GetVoltageCheck() && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );

		return true;
	}

    double WaveguideModes::GetVoltagePhaseAntenna(int aChannel, int l, int m, int n)
    {
    	return fVoltagePhaseAntenna[aChannel][l][m][n];
    }
    void WaveguideModes::SetVoltagePhaseAntenna( double aPhase, int aChannel, int l, int m, int n )
    {
        fVoltagePhaseAntenna[aChannel][l][m][n] = aPhase;
    }
    double WaveguideModes::GetVoltagePhaseShort(int aChannel, int l, int m, int n)
    {
    	return fVoltagePhaseShort[aChannel][l][m][n];
    }
    void WaveguideModes::SetVoltagePhaseShort( double aPhase, int aChannel, int l, int m, int n )
    {
        fVoltagePhaseShort[aChannel][l][m][n] = aPhase;
    }




} /* namespace locust */

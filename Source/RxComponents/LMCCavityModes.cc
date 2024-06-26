/*
 * LMCCavityModes.cc
 *
 *  Created on: Jul 16, 2021
 *      Author: pslocum
 */

#include "LMCCavityModes.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "CavityModes" );

    CavityModes::CavityModes():
		fVoltagePhase( 0 ),
		fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    CavityModes::~CavityModes()
    {
    }


    bool CavityModes::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombiner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from CavityModes subclass");
    		return false;
    	}

        SetNCavityModes(fInterface->fField->GetNModes());

        SizeNChannels(GetNChannels());

    	return true;
    }

    bool CavityModes::SizeNChannels(int aNumberOfChannels)
    {
    	SetNChannels(aNumberOfChannels);

    	std::vector<std::vector<std::vector<std::vector<double>>>> tZeroVector;
    	fVoltagePhase.swap(tZeroVector);
    	fVoltagePhase.resize(aNumberOfChannels);

    	for (int n = 0; n < GetNChannels(); n++)
    	{
    		fVoltagePhase[n].resize(GetNCavityModes());
    		for (int i = 0; i < GetNCavityModes(); i++)
    		{
    			fVoltagePhase[n][i].resize(GetNCavityModes());
    			for (int j = 0; j < GetNCavityModes(); j++)
    			{
    				fVoltagePhase[n][i][j].resize(GetNCavityModes());
    			}
    		}
    	}

    	return true;
    }

	bool CavityModes::AddOneModeToCavityProbe(int l, int m, int n, Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> cavityDopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, int channelIndex, bool initParticle)
	{
		double dopplerFrequency = cavityDopplerFrequency[0];  // Only one shift, unlike in waveguide.
		SetVoltagePhase( GetVoltagePhase(channelIndex, l, m, n) + dopplerFrequency * dt, channelIndex, l, m, n ) ;
		double voltageValue = excitationAmplitude * EFieldAtProbe;
		voltageValue *= cos(GetVoltagePhase(channelIndex, l, m, n));

		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);

		if ( GetVoltageCheck() && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );
		return true;
	}

    double CavityModes::GetVoltagePhase(int aChannel, int l, int m, int n)
    {
    	return fVoltagePhase[aChannel][l][m][n];
    }

    void CavityModes::SetVoltagePhase ( double aPhase, int aChannel, int l, int m, int n )
    {
        fVoltagePhase[aChannel][l][m][n] = aPhase;
    }


} /* namespace locust */

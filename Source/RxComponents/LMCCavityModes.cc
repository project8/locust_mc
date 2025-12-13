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
        fVoltagePhase( {{{{0.0}}}} ),
        fModeMap( false ),
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

        if( aParam.has( "upload-modemap-filename" ) )
        {
            fModeMap = true;
        }

        fModeSet = fInterface->fField->ModeSelect(fInterface->fbWaveguide, 0);

    	return true;
    }

    bool CavityModes::SetChannelPhaseShifts()
    {
        if (( GetNChannels() > 3 ) || ( fModeSet.size() != GetNChannels() ) || GetNChannels() != fInterface->fField->GetNChannels())
        {
            LERROR(lmclog,"LMCCavityModes sees " << GetNChannels() << " channels and " << fModeSet.size() << " modes.");
            LERROR(lmclog,"The cavity simulation only supports up to 3 channels right now, and the "
    	    		"number of channels has to be the same as the number of modes.");
            exit(-1);
        }

        fChannelPhaseShifts.resize(GetNChannels());
        for (int channelIndex = 0; channelIndex < GetNChannels(); channelIndex++)
        {
            fChannelPhaseShifts[channelIndex].resize(fModeSet.size());
            for (int mu = 0; mu < fModeSet.size(); mu++)
            {
                int modeSignThetaComp = 1;
                double signalPhaseShift = 0.;
                if ( !fModeMap )
                {
                    bool bTE = fModeSet[mu][0];
                    int l = fModeSet[mu][1];
                    int m = fModeSet[mu][2];
                    int n = fModeSet[mu][3];

                    std::vector<double> emptyPositionVector = {0., 0., 0.};
                    modeSignThetaComp = ( fInterface->fField->GetFieldAtProbe(l,m,n,0,emptyPositionVector,1)[channelIndex] < 0. );
                    // phase shift of PI for negative mode field at probe:
                    fChannelPhaseShifts[channelIndex][mu] = LMCConst::Pi() * ( modeSignThetaComp );
                    LPROG(lmclog,"Channel "<< channelIndex << " phase shift for mode " <<bTE<<l<<m<<n<< " is " << fChannelPhaseShifts[channelIndex][mu]);
                }
            }
        }

        LPROG( lmclog, "Channel phase offsets have been calculated." );

        return true;
    }


    bool CavityModes::SizeNChannels(int aNumberOfChannels)
    {
        SetNChannels(aNumberOfChannels);
        SetChannelPhaseShifts();

        std::vector<std::vector<double>> tZeroVector = {{0.}};
        fVoltagePhase.swap(tZeroVector);
        fVoltagePhase.resize(aNumberOfChannels);

        for (int channelIndex = 0; channelIndex < GetNChannels(); channelIndex++)
        {
            fVoltagePhase[channelIndex].resize(fModeSet.size());
            for (int mu = 0; mu < fModeSet.size(); mu++)
            {
    	        fVoltagePhase[channelIndex][mu] = fChannelPhaseShifts[channelIndex][mu];
            }
        }

        LPROG( lmclog, "Channel phase offsets have been applied." );

        return true;
    }

	bool CavityModes::AddOneModeToCavityProbe(int mu, Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> cavityDopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, int channelIndex, bool initParticle)
	{

		double dopplerFrequency = cavityDopplerFrequency[0];  // Only one shift, unlike in waveguide.
		SetVoltagePhase( GetVoltagePhase(channelIndex, mu) + dopplerFrequency * dt, channelIndex, mu ) ;
		double voltageValue = excitationAmplitude * EFieldAtProbe;
		voltageValue *= cos(GetVoltagePhase(channelIndex, mu) + fChannelPhaseShifts[channelIndex][mu] );

		aSignal->LongSignalTimeComplex()[sampleIndex][0] += 2. * voltageValue * totalScalingFactor * sin(phi_LO);
		aSignal->LongSignalTimeComplex()[sampleIndex][1] += 2. * voltageValue * totalScalingFactor * cos(phi_LO);

		if ( GetVoltageCheck() && (sampleIndex%100 < 1) )
			LPROG( lmclog, "Voltage " << sampleIndex << " is <" << aSignal->LongSignalTimeComplex()[sampleIndex][1] << ">" );
		return true;
	}

    double CavityModes::GetVoltagePhase(int aChannel, int mu)
    {
    	return fVoltagePhase[aChannel][mu];
    }

    void CavityModes::SetVoltagePhase ( double aPhase, int aChannel, int mu )
    {
        fVoltagePhase[aChannel][mu] = aPhase;
    }

} /* namespace locust */

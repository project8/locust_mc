/*
 * LMCAntennaSignalTransmitter.cc
 *
 *  Created on: May 4, 2019
 *      Author: pslocum
 */

#include "LMCAntennaSignalTransmitter.hh"

#include "logger.hh"
using std::string;


namespace locust
{
    LOGGER( lmclog, "AntennaSignalTransmitter" );
    
    AntennaSignalTransmitter::AntennaSignalTransmitter() :
    fInputSignalType(1),
    fInputFrequency( 0.0 ),
    fInputAmplitude(1.0),
    fAntennaType(0)
    {
    }
    
    AntennaSignalTransmitter::~AntennaSignalTransmitter()
    {
    }
    
    bool AntennaSignalTransmitter::Configure( const scarab::param_node& aParam )
    {
        if(!fTransmitterHandler.Configure(aParam))
        {
            LERROR(lmclog,"Error configuring FIRHandler class");
        }
        
        if( aParam.has( "input-signal-type" ) )
        {
            fInputSignalType = aParam["input-signal-type"]().as_int();
        }
        
	if( aParam.has( "transmitter-frequency" ) )
        {
            fInputFrequency= aParam["transmitter-frequency"]().as_double();
        }
        
        if( aParam.has( "antenna-voltage-amplitude" ) )
        {
            fInputAmplitude = aParam["antenna-voltage-amplitude"]().as_double();
        }

        if( aParam.has( "transmitter-antenna-type" ) )
        {
            fAntennaType = SetAntennaType(aParam["transmitter-antenna-type"]().as_string());
        }
        else
        {
            fAntennaType = SetAntennaType("antenna-signal-dipole");
        }

     	if(!fTransmitterHardware->Configure(aParam))
	    {
     		LERROR(lmclog,"Error configuring TransmitterHardware class");
	    }


        return true;
    }
    
    bool AntennaSignalTransmitter::SetAntennaType( std::string transmitterAntennaType )
     {
     	if (transmitterAntennaType == "antenna-signal-dipole")
     		{
     			fAntennaType = 0; // default
     			fTransmitterHardware = new DipoleAntenna;
     		}
     	else if (transmitterAntennaType == "antenna-signal-turnstile")
     		{
     			fAntennaType = 1;
     			fTransmitterHardware = new TurnstileAntenna;
     		}
     	else return false;

     	return true;
     }


    double* AntennaSignalTransmitter::GetEFieldCoPol(int fieldPointIndex, double dt)
    {
    	LMCThreeVector pointOfInterest=GetFieldPoint(fieldPointIndex);
        double estimatedField=0.0;
        //if ( ( zIndex == 0 ) && (channelIndex == 0) ) fPhaseDelay+= 2.*LMCConst::Pi()*fInputFrequency*dt;
        if (fieldPointIndex == 0) fPhaseDelay+= 2.*LMCConst::Pi()*fInputFrequency*dt;
        double voltagePhase=fPhaseDelay - GetPropagationPhaseDelay(fieldPointIndex);

        for (unsigned iAntenna = 0; iAntenna < fTransmitterHardware->GetNAntennas(); iAntenna++)
        {
        	// phase shift between transmitting antennas if nAntennas > 1.
        	voltagePhase += iAntenna * fTransmitterHardware->GetDrivePhaseDifference();

        	if(fInputSignalType==1) //sinusoidal wave for dipole antenna
        	{
        		for( unsigned index = 0; index <fTransmitterHandler.GetFilterSize();index++)
        		{
        			double voltageValue = GetFieldAtOrigin(fInputAmplitude,voltagePhase);
        			delayedVoltageBuffer[0].push_back(voltageValue);
        			delayedVoltageBuffer[0].pop_front();

        			voltagePhase += 2.*LMCConst::Pi()*fInputFrequency*fTransmitterHandler.GetFilterResolution();
        		}
        	}
        
        	else// For now using sinusoidal as well
        	{
        		for( unsigned index = 0; index <fTransmitterHandler.GetFilterSize();index++)
        		{
        			double voltageValue = GetFieldAtOrigin(fInputAmplitude,voltagePhase);
        			delayedVoltageBuffer[0].push_back(voltageValue);
        			delayedVoltageBuffer[0].pop_front();
        			voltagePhase += 2.*LMCConst::Pi()*fInputFrequency*fTransmitterHandler.GetFilterResolution();
        		}
        	}

        	// find total field from all transmitting antennas in fTransmitterHardware object.
        	estimatedField += fTransmitterHandler.ConvolveWithFIRFilter(delayedVoltageBuffer[0]) * fTransmitterHardware->GetPatternFactor(pointOfInterest, iAntenna);

        } // nAntennas

        double* FieldSolution = new double[2];
        FieldSolution[0] = estimatedField / fTransmitterHardware->GetPropagationDistance(pointOfInterest); // field at point
        FieldSolution[1] = 2. * LMCConst::Pi() * fInputFrequency; // rad/s

        return FieldSolution;
    }
    
    bool AntennaSignalTransmitter::InitializeTransmitter()
    {
        if(!fTransmitterHandler.ReadHFSSFile())
        {
	    LERROR(lmclog,"Error reading HFSS file");
            return false;
        }
        double filterSize=fTransmitterHandler.GetFilterSize();
        InitializeBuffers(filterSize);
        fInitialPhaseDelay = -2.*LMCConst::Pi()*(filterSize*fTransmitterHandler.GetFilterResolution())*fInputFrequency;
        fPhaseDelay = fInitialPhaseDelay;
        return true;
    }
    
    void AntennaSignalTransmitter::AddIncidentKVector(LMCThreeVector pointOfInterest)
    {
	LMCThreeVector incidentKVector=fTransmitterHardware->ExtractIncidentKVector(pointOfInterest); 
	Transmitter::AddIncidentKVector(incidentKVector);
    }
 
    double AntennaSignalTransmitter::GetInitialPhaseDelay()
    {
        return fInitialPhaseDelay;
    }
    
    void AntennaSignalTransmitter::AddPropagationPhaseDelay(LMCThreeVector pointOfInterest)
    {
        double phaseDelay = 2.*LMCConst::Pi()*fInputFrequency/LMCConst::C()*
		fTransmitterHardware->GetPropagationDistance(pointOfInterest);
	Transmitter::AddPropagationPhaseDelay(phaseDelay);
    }

    void AntennaSignalTransmitter::InitializeBuffers(unsigned filterbuffersize)
    {
        FieldBuffer aFieldBuffer;
        delayedVoltageBuffer = aFieldBuffer.InitializeBuffer(1,1,filterbuffersize);
    }
    
    double AntennaSignalTransmitter::GetFieldAtOrigin(double inputAmplitude,double voltagePhase)
    {
        //double normalizedVoltage = cos(voltagePhase);
        double normalizedDerivative = ApplyDerivative(voltagePhase);
        // Only missing tau, f_g
        // And distance will be applied later
        double field = inputAmplitude*normalizedDerivative/(2*LMCConst::Pi()*LMCConst::C());
        return field;
    }
} /* namespace locust */

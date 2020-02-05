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
    fAntennaPositionX( 0.0 ),
    fAntennaPositionY( 0.0 ),
    fAntennaPositionZ( 0.0 ),
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

        if( aParam.has( "antenna-x-position" ) )
        {
            fAntennaPositionX= aParam["antenna-x-position"]().as_double();
        }
        
        if( aParam.has( "antenna-y-position" ) )
        {
            fAntennaPositionY = aParam["antenna-y-position"]().as_double();
        }
        
        if( aParam.has( "antenna-z-position" ) )
        {
            fAntennaPositionZ = aParam["antenna-z-position"]().as_double();
        }
        
        if( aParam.has( "antenna-voltage-amplitude" ) )
        {
            fInputAmplitude = aParam["antenna-voltage-amplitude"]().as_double();
        }

        if( aParam.has( "transmitter-antenna-type" ) )
        {
            fAntennaType = SetAntennaType(aParam["transmitter-antenna-type"]().as_string());
        }

        return true;
    }
    
    bool AntennaSignalTransmitter::SetAntennaType( std::string transmitterAntennaType )
     {
     	if (transmitterAntennaType == "antenna-signal-dipole") fAntennaType = 0; // default
     	else if (transmitterAntennaType == "antenna-signal-turnstile") fAntennaType = 1;
     	else return false;
     	return true;
     }


    LMCThreeVector AntennaSignalTransmitter::GetAntennaPosition() const
    {
        return fAntennaPosition;
    }
    
    void AntennaSignalTransmitter::SetAntennaPosition(const LMCThreeVector &antennaPosition)
    {
        fAntennaPosition=antennaPosition;
    }


    double AntennaSignalTransmitter::GetPropagationDistance(Receiver* currentElement)
    {
        double relativePatchPosX=currentElement->GetPosition().GetX() - fAntennaPosition.GetX();
        double relativePatchPosY=currentElement->GetPosition().GetY() - fAntennaPosition.GetY();
        double relativePatchPosZ=currentElement->GetPosition().GetZ() - fAntennaPosition.GetZ();
        double propagationDistance = sqrt(relativePatchPosX*relativePatchPosX+relativePatchPosY*relativePatchPosY+relativePatchPosZ*relativePatchPosZ);
        return propagationDistance;
    }


    double AntennaSignalTransmitter::GetPropagationPhaseChange(Receiver* currentElement)
    {
        double relativePatchPosX=currentElement->GetPosition().GetX() - fAntennaPosition.GetX();
        double relativePatchPosY=currentElement->GetPosition().GetY() - fAntennaPosition.GetY();
        double relativePatchPosZ=currentElement->GetPosition().GetZ() - fAntennaPosition.GetZ();
        double phaseChange = 2.*LMCConst::Pi()*fInputFrequency/LMCConst::C()*GetPropagationDistance(currentElement);
    	return phaseChange;
    }

    double AntennaSignalTransmitter::GetAOIFactor(Receiver* currentElement)
    {
        double relativePatchPosX=currentElement->GetPosition().GetX() - fAntennaPosition.GetX();
        double relativePatchPosY=currentElement->GetPosition().GetY() - fAntennaPosition.GetY();
        double relativePatchPosZ=currentElement->GetPosition().GetZ() - fAntennaPosition.GetZ();
        double elementAntennaDistance = (currentElement->GetPosition() - fAntennaPosition).Magnitude();

        double aoiFactor = (relativePatchPosX*currentElement->GetNormalDirection().GetX() +
        		relativePatchPosY*currentElement->GetNormalDirection().GetY() +
				relativePatchPosZ*currentElement->GetNormalDirection().GetZ()) /
						elementAntennaDistance;

    	return aoiFactor;
    }

    double* AntennaSignalTransmitter::GetEFieldCoPol(Receiver* currentElement, int channelIndex, int zIndex, double elementSpacing, int nElementsPerStrip, double dt)
    {
        double estimatedField=0.0;
        if ( ( zIndex == 0 ) && (channelIndex == 0) ) fPhaseDelay+= 2.*LMCConst::Pi()*fInputFrequency*dt;
        double voltagePhase=fPhaseDelay + GetPropagationPhaseChange(currentElement);

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

        estimatedField=fTransmitterHandler.ConvolveWithFIRFilter(delayedVoltageBuffer[0]) * GetAOIFactor(currentElement);
        double* FieldSolution = new double[2];
        FieldSolution[0] = estimatedField / GetPropagationDistance(currentElement); // field at Rx antenna.
        FieldSolution[1] = 2. * LMCConst::Pi() * fInputFrequency; // rad/s

        return FieldSolution;
    }
    
    bool AntennaSignalTransmitter::InitializeTransmitter()
    {
        fAntennaPosition.SetComponents(fAntennaPositionX,fAntennaPositionY,fAntennaPositionZ);
        
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
    
    double AntennaSignalTransmitter::GetInitialPhaseDelay()
    {
        return fInitialPhaseDelay;
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

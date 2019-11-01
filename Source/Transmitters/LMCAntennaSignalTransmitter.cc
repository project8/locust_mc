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
    fAntennaPositionX( 0.0 ),
    fAntennaPositionY( 0.0 ),
    fAntennaPositionZ( 0.0 ),
    fInputAmplitude(1)
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
        
        if( aParam.has( "input-signal-frequency" ) )
        {
            fInputFrequency= aParam["input-signal-frequency"]().as_double();
        }
        
        if( aParam.has( "array-radius" ) )
        {
            fArrayRadius = aParam["array-radius"]().as_double();
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
        
        if( aParam.has( "input-signal-amplitude" ) )
        {
            fInputAmplitude = aParam["input-signal-amplitude"]().as_double();
        }
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
    double AntennaSignalTransmitter::GenerateSignal(Signal *aSignal,double acquisitionRate)
    {
        double estimatedField=0.0;
        double voltagePhase=fPhaseDelay;
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
        estimatedField=fTransmitterHandler.ConvolveWithFIRFilter(delayedVoltageBuffer[0]);
        fPhaseDelay+= 2.*LMCConst::Pi()*fInputFrequency/aSignal->DecimationFactor()/(acquisitionRate*1.e6);
        return estimatedField;
    }
    
    bool AntennaSignalTransmitter::InitializeTransmitter()
    {
        fAntennaPosition.SetComponents(fAntennaPositionX,fAntennaPositionY,fAntennaPositionZ);
        
        if(!fTransmitterHandler.ReadHFSSFile())
        {
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

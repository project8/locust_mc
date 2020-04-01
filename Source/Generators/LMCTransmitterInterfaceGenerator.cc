/*
 * LMCTransmitterInterfaceGenerator.cc
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
 */

#include "LMCTransmitterInterfaceGenerator.hh"
#include "LMCEventHold.hh"
#include "LMCRunKassiopeia.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>

#include "LMCDigitizer.hh"
#include <chrono>


namespace locust
{
    LOGGER( lmclog, "TransmitterInterfaceGenerator" );

    MT_REGISTER_GENERATOR(TransmitterInterfaceGenerator, "transmitter-interface");

    TransmitterInterfaceGenerator::TransmitterInterfaceGenerator( const std::string& aName ) :
        Generator( aName ),
        gxml_filename("blank.xml"),
	fTextFileWriting( 0 ),
        EFieldBuffer( 1 ),
        EFrequencyBuffer( 1 ),
        fFieldBufferSize( 50 ),
        fdtFilter( 50 ),
	fSwapFrequency( 1000 ),
	fNPoints(50)
    {
        fRequiredSignalState = Signal::kTime;
    }

    TransmitterInterfaceGenerator::~TransmitterInterfaceGenerator()
    {
    }

    bool TransmitterInterfaceGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "transmitter" ))
        {
        	int ntransmitters = 0;
        	if(aParam["transmitter"]().as_string() == "antenna")
        	{
        		ntransmitters += 1;
			fTransmitter = std::make_shared<AntennaSignalTransmitter>();
        		if(!fTransmitter->Configure(aParam))
        		{
        			LERROR(lmclog,"Error Configuring antenna signal transmitter class");
        		}
        		if(!fTransmitter->InitializeTransmitter())
        		{
        			exit(-1);
        		}
        	}

        	if(aParam["transmitter"]().as_string() == "planewave")
        	{
        		ntransmitters += 1;
			fTransmitter = std::make_shared<PlaneWaveTransmitter>();
        		if(!fTransmitter->Configure(aParam))
        		{
        			LERROR(lmclog,"Error Configuring planewave transmitter class");
        		}
        	}

        	if(aParam["transmitter"]().as_string() == "kassiopeia")
        	{
        		ntransmitters += 1;
			fTransmitter = std::make_shared<KassTransmitter>();
        		if(!fTransmitter->Configure(aParam))
        		{
        			LERROR(lmclog,"Error Configuring kassiopeia transmitter class");
        		}
        	}

        	if (ntransmitters != 1)
        	{
        		LERROR(lmclog,"The generator needs a single transmitter.  Please choose transmitter:antenna or transmitter:planewave or transmitter:kassiopeia in the config file.");
                exit(-1);
        	}
        }
        else
        {
    		LERROR(lmclog,"The generator has been configured without a transmitter.  Please choose transmitter:antenna or transmitter:planewave or transmitter:kassiopeia in the config file.");
            exit(-1);
        }

        if( aParam.has( "buffer-size" ) )
        {
        	fFieldBufferSize = aParam["buffer-size"]().as_int();
        }

        if( aParam.has( "swap-frequency" ) )
        {
            fSwapFrequency = aParam["swap-frequency"]().as_int();
        }

        if( aParam.has( "xml-filename" ) )
        {
            gxml_filename = aParam["xml-filename"]().as_string();
        }
        if( aParam.has( "text-filewriting" ) )
        {
            fTextFileWriting = aParam["text-filewriting"]().as_bool();
        }
        return true;
    }

    void TransmitterInterfaceGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    void TransmitterInterfaceGenerator::InitializeFieldPoint(LMCThreeVector pointOfInterest)
    {
            fTransmitter->InitializeFieldPoint(pointOfInterest);
	    fAllFieldCopol.push_back(pointOfInterest);
    }

    bool TransmitterInterfaceGenerator::WakeBeforeEvent()
    {
        fPreEventCondition.notify_one();
        return true;
    }

    bool TransmitterInterfaceGenerator::ReceivedKassReady()
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        printf("LMC about to wait ..\n");

        if( !fKassEventReady)
        {
            std::unique_lock< std::mutex >tLock( fKassReadyMutex );
            fKassReadyCondition.wait( tLock );
        }

        if (fFalseStartKassiopeia)  // workaround for some Macs
        {
            std::unique_lock< std::mutex >tLock( fKassReadyMutex );
            fKassReadyCondition.wait( tLock );
        }

        return true;
    }

    // fields incident on element.
    void TransmitterInterfaceGenerator::RecordIncidentFields(FILE *fp,  double t_old,LMCThreeVector pointOfInterest, double tEFieldCoPol)
    {
    	if (t_old > 0.5e-9)
         	{
         	fprintf(fp, "%d %g %g\n", pointOfInterest.GetX(), pointOfInterest.GetY(),pointOfInterest.GetZ(),tEFieldCoPol);
         	}
    }

    void TransmitterInterfaceGenerator::FillBuffers(Signal* aSignal, double DopplerFrequency, double EFieldValue,unsigned pointIndex, unsigned timeIndex)
    {
    	EFieldBuffer[pointIndex].push_back(EFieldValue);
    	EFrequencyBuffer[pointIndex].push_back(DopplerFrequency/2./LMCConst::Pi());
    }

    void TransmitterInterfaceGenerator::PopBuffers(unsigned pointIndex)
    {
    	EFieldBuffer[pointIndex].pop_front();
    	EFrequencyBuffer[pointIndex].pop_front();
    }

    void TransmitterInterfaceGenerator::InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize)
    {
    	EFieldBuffer = fFieldBuffer.InitializeBuffer(fNPoints, fieldbuffersize);
    	EFrequencyBuffer = fFieldBuffer.InitializeBuffer(fNPoints, fieldbuffersize);
    }

    void TransmitterInterfaceGenerator::CleanupBuffers()
    {
    	EFieldBuffer = fFieldBuffer.CleanupBuffer(EFieldBuffer);
    	EFrequencyBuffer = fFieldBuffer.CleanupBuffer(EFrequencyBuffer);
    }

} /* namespace locust */


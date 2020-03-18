/*
 * LMCFreeSpaceGenerator.cc
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
 */

#include "LMCFreeSpaceGenerator.hh"
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
    LOGGER( lmclog, "FreeSpaceGenerator" );

    MT_REGISTER_GENERATOR(FreeSpaceGenerator, "free-space");

    FreeSpaceGenerator::FreeSpaceGenerator( const std::string& aName ) :
        Generator( aName ),
        gxml_filename("blank.xml"),
	fTextFileWriting( 0 ),
        EFieldBuffer( 1 ),
        EFrequencyBuffer( 1 ),
        IndexBuffer( 1 ),
        fFieldBufferSize( 50 ),
        fNFilterBins(50),
        fdtFilter( 50 ),
	fSwapFrequency( 1000 )
    {
	std::cout<< "FreeSpaceGenerator::FreeSpaceGenerator "<<std::endl;
        fRequiredSignalState = Signal::kTime;
    }

    FreeSpaceGenerator::~FreeSpaceGenerator()
    {
	std::cout<< "FreeSpaceGenerator::~FreeSpaceGenerator "<<std::endl;
    }

    bool FreeSpaceGenerator::Configure( const scarab::param_node& aParam )
    {
	std::cout<< "FreeSpaceGenerator::Configure"<<std::endl;
        if( aParam.has( "transmitter" ))
        {
        	int ntransmitters = 0;
        	if(aParam["transmitter"]().as_string() == "antenna")
        	{
        		ntransmitters += 1;
			fTransmitter = new AntennaSignalTransmitter;
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
			fTransmitter = new PlaneWaveTransmitter;
        		if(!fTransmitter->Configure(aParam))
        		{
        			LERROR(lmclog,"Error Configuring planewave transmitter class");
        		}
        	}

        	if(aParam["transmitter"]().as_string() == "kassiopeia")
        	{
        		ntransmitters += 1;
			fTransmitter = new KassTransmitter;
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

    void FreeSpaceGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
	std::cout<< "FreeSpaceGenerator::Accept"<<std::endl;
        aVisitor->Visit( this );
        return;
    }


    static void* KassiopeiaInit(const std::string &aFile)
    {
	std::cout<< "FreeSpaceGenerator::KassiopeiaInit"<<std::endl;
        RunKassiopeia RunKassiopeia1;
        RunKassiopeia1.Run(aFile);
        RunKassiopeia1.~RunKassiopeia();

        return 0;
    }

    void FreeSpaceGenerator::InitializeFieldPoints()
    {
	std::cout<< "FreeSpaceGenerator::InitializeFieldPoints"<<std::endl;
	fNPoints=50;
	for(int pointIndex = 0; pointIndex< fNPoints; ++pointIndex)
	{
	    LMCThreeVector point(0.0,0.0,0.0);
            fTransmitter->InitializeFieldPoint(point);
	    fAllFieldCopol[pointIndex]=point;
	}
    }

    bool FreeSpaceGenerator::WakeBeforeEvent()
    {
	std::cout<< "FreeSpaceGenerator::WakeBeforeEvent"<<std::endl;
        fPreEventCondition.notify_one();
        return true;
    }

    bool FreeSpaceGenerator::ReceivedKassReady()
    {
	std::cout<< "FreeSpaceGenerator::ReceivedKassReady"<<std::endl;
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
    void FreeSpaceGenerator::RecordIncidentFields(FILE *fp,  double t_old, LMCThreeVector pointOfInterest, double tEFieldCoPol)
    {
	std::cout<< "FreeSpaceGenerator::RecordIncidentFields"<<std::endl;
    	if (t_old > 0.5e-9)
         	{
         	fprintf(fp, "%d %g %g\n", pointOfInterest.GetX(), pointOfInterest.GetY(),pointOfInterest.GetZ(),tEFieldCoPol);
         	}
	std::cout<< pointOfInterest.GetX() << " : "<< pointOfInterest.GetY() << " : "<< pointOfInterest.GetZ()<< " : "<< tEFieldCoPol<<std::endl;
    }

    void FreeSpaceGenerator::DriveAntenna(FILE *fp, int PreEventCounter, unsigned index, Signal* aSignal, int nfilterbins, double dtfilter)
    {
	std::cout<< "FreeSpaceGenerator::DriveAntenna"<<std::endl;
        const int signalSize = aSignal->TimeSize();
        unsigned pointIndex = 0;
        unsigned sampleIndex = 0;

        unsigned tTotalElementIndex = 0;

        for(int pointIndex = 0; pointIndex < fNPoints; ++pointIndex)
        {
                sampleIndex = pointIndex*aSignal->DecimationFactor() + index;  // which point and which sample

                double* tFieldSolution = new double[2];
                if (!fTransmitter->IsKassiopeia())
                {
                	tFieldSolution = fTransmitter->GetEFieldCoPol(pointIndex, 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor()));
                }
                else
                {
                	tFieldSolution = fTransmitter->SolveKassFields(fAllFieldPoints[pointIndex],fAllFieldCopol[pointIndex],t_old,pointIndex);
                }
                if (fTextFileWriting==1) {}
		RecordIncidentFields(fp, t_old,fAllFieldPoints[pointIndex], tFieldSolution[1]);
 	        FillBuffers(aSignal, tFieldSolution[1], tFieldSolution[0],pointIndex,index);
                PopBuffers(pointIndex);
        } // channels loop

        t_old += 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        if ( index%fSwapFrequency == 0 ) CleanupBuffers();  // release memory
    }

    void FreeSpaceGenerator::FillBuffers(Signal* aSignal, double DopplerFrequency, double EFieldValue,unsigned pointIndex, unsigned timeIndex)
    {
    	EFieldBuffer[pointIndex].push_back(EFieldValue);
    	EFrequencyBuffer[pointIndex].push_back(DopplerFrequency/2./LMCConst::Pi());
    	IndexBuffer[pointIndex].push_back(pointIndex*aSignal->TimeSize()*aSignal->DecimationFactor() + timeIndex);
    }

    void FreeSpaceGenerator::PopBuffers(unsigned pointIndex)
    {
    	EFieldBuffer[pointIndex].pop_front();
    	EFrequencyBuffer[pointIndex].pop_front();
    	IndexBuffer[pointIndex].pop_front();
    }

    void FreeSpaceGenerator::InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize)
    {
    	FieldBuffer aFieldBuffer;
    	EFieldBuffer = aFieldBuffer.InitializeBuffer(fNPoints, fieldbuffersize);
    	EFrequencyBuffer = aFieldBuffer.InitializeBuffer(fNPoints, fieldbuffersize);
    	IndexBuffer = aFieldBuffer.InitializeUnsignedBuffer(fNPoints, fieldbuffersize);
    }

    void FreeSpaceGenerator::CleanupBuffers()
    {
    	FieldBuffer aFieldBuffer;
    	EFieldBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    	EFrequencyBuffer = aFieldBuffer.CleanupBuffer(EFrequencyBuffer);
    	IndexBuffer = aFieldBuffer.CleanupBuffer(IndexBuffer);
    }

    bool FreeSpaceGenerator::DoGenerate( Signal* aSignal )
    {
	std::cout<< "FreeSpaceGenerator::DriveAntenna"<<std::endl;
        FILE *fp = fopen("incidentfields.txt", "w");

        //n samples for event spacing in Kass.
        int PreEventCounter = 0;
        const int NPreEventSamples = 150000;
        fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        InitializeBuffers(fNFilterBins, fFieldBufferSize);
        InitializeFieldPoints();

        if (!fTransmitter->IsKassiopeia())
        {
        	for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        	{
        		DriveAntenna(fp, PreEventCounter, index, aSignal, fNFilterBins,fdtFilter);
        	}  // for loop
        	return true;
        }

        else if (fTransmitter->IsKassiopeia())
        {

            std::thread Kassiopeia(KassiopeiaInit, gxml_filename);  // spawn new thread
        	fRunInProgress = true;

        for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        {
            if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
            {
                if (ReceivedKassReady()) fPreEventInProgress = true;
            }

            if (fPreEventInProgress)
            {
                PreEventCounter += 1;
                //printf("preeventcounter is %d\n", PreEventCounter);
                if (PreEventCounter > NPreEventSamples)  // finished noise samples.  Start event.
                {
                    fPreEventInProgress = false;  // reset.
                    fEventInProgress = true;
                    //printf("LMC about to wakebeforeevent\n");
                    WakeBeforeEvent();  // trigger Kass event.
                }
            }

            if (fEventInProgress)  // fEventInProgress
                if (fEventInProgress)  // check again.
                {
                    //printf("waiting for digitizer trigger ... index is %d\n", index);
                    std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );
                    tLock.lock();
                    fDigitizerCondition.wait( tLock );
                    if (fEventInProgress)
                    {
                        //printf("about to drive antenna, PEV is %d\n", PreEventCounter);
                        DriveAntenna(fp, PreEventCounter, index, aSignal, fNFilterBins,fdtFilter);
                        PreEventCounter = 0; // reset
                    }
                    tLock.unlock();
                }

        }  // for loop

        printf("finished signal loop\n");

        fclose(fp);
        fRunInProgress = false;  // tell Kassiopeia to finish.
        fDoneWithSignalGeneration = true;  // tell LMCCyclotronRadExtractor
        WakeBeforeEvent();
        Kassiopeia.join();

        return true;
        }

    }

} /* namespace locust */


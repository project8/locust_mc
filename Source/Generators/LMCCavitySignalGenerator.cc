/*
 * LMCCavitySignalGenerator.cc
 *
 *  Created on: Mar 30, 2021
 *      Author: pslocum
 */

#include "LMCCavitySignalGenerator.hh"

#include "LMCRunKassiopeia.hh"

#include "logger.hh"

#include <chrono>
#include <thread>


namespace locust
{
    LOGGER( lmclog, "CavitySignalGenerator" );

    MT_REGISTER_GENERATOR(CavitySignalGenerator, "cavity-signal");

    CavitySignalGenerator::CavitySignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fLO_Frequency( 0.),
        fArrayRadius( 0. ),
        fNElementsPerStrip( 0. ),
		fNSubarrays( 1 ),
		fZShiftArray( 0. ),
        fElementSpacing( 0. ),
        gxml_filename("blank.xml"),
		fTextFileWriting( 0 ),
        fphiLO(0.),
		fNPreEventSamples( 150000 ),
		fThreadCheckTime(200),
		fKassNeverStarted( false ),
		fSkippedSamples( false ),
		fInterface( new KassLocustInterface() )
    {
        fRequiredSignalState = Signal::kTime;

        KLInterfaceBootstrapper::get_instance()->SetInterface( fInterface );
    }

    CavitySignalGenerator::~CavitySignalGenerator()
    {
    }

    bool CavitySignalGenerator::Configure( const scarab::param_node& aParam )
    {

        if( aParam.has( "transmitter" ))
        {
        	int ntransmitters = 0;

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
        		LERROR(lmclog,"LMCCavitySignalGenerator needs a single transmitter.  Please choose transmitter:antenna or transmitter:planewave or transmitter:kassiopeia in the config file.");
                exit(-1);
        	}
        }
        else
        {
    		LERROR(lmclog,"LMCCavitySignalGenerator has been configured without a transmitter.  Please choose transmitter:antenna or transmitter:planewave or transmitter:kassiopeia in the config file.");
            exit(-1);
        }

        if( aParam.has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam["lo-frequency"]().as_double();
        }

        if( aParam.has( "array-radius" ) )
        {
            fArrayRadius = aParam["array-radius"]().as_double();
        }

        if( aParam.has( "nelements-per-strip" ) )
        {
            fNElementsPerStrip = aParam["nelements-per-strip"]().as_int();
        }

        if( aParam.has( "n-subarrays" ) )
        {
            fNSubarrays = aParam["n-subarrays"]().as_int();
        }

        if( aParam.has( "element-spacing" ) )
        {
            fElementSpacing = aParam["element-spacing"]().as_double();
        }
        if( aParam.has( "zshift-array" ) )
        {
            fZShiftArray = aParam["zshift-array"]().as_double();
        }
        if( aParam.has( "event-spacing-samples" ) )
        {
            fNPreEventSamples = aParam["event-spacing-samples"]().as_int();
        }
        if( aParam.has( "thread-check-time" ) )
        {
            fThreadCheckTime = aParam["thread-check-time"]().as_int();
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

    void CavitySignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }



    void CavitySignalGenerator::KassiopeiaInit(const std::string &aFile)
    {
        RunKassiopeia tRunKassiopeia;
        tRunKassiopeia.Run(aFile, fInterface);
        return;
    }

    void CavitySignalGenerator::InitializeFieldPoints(std::vector< Channel<Receiver*> > allRxChannels)
    {
	for(int channelIndex = 0; channelIndex < fNChannels; ++channelIndex)
	{
            for(int elementIndex = 0; elementIndex < fNElementsPerStrip; ++elementIndex)
            {
            	fTransmitter->InitializeFieldPoint(allRxChannels[channelIndex][elementIndex]->GetPosition());
            }
	}
    }

    void CavitySignalGenerator::WakeBeforeEvent()
    {
        fInterface->fPreEventCondition.notify_one();
        return;
    }

    bool CavitySignalGenerator::ReceivedKassReady()
    {

    	std::this_thread::sleep_for(std::chrono::milliseconds(100));
        LPROG( lmclog, "LMC about to wait" );

        if(!fInterface->fKassEventReady)
        {
            std::unique_lock< std::mutex >tLock( fInterface->fKassReadyMutex );
            fInterface->fKassReadyCondition.wait( tLock );
            return true;
        }
        else if (fInterface->fKassEventReady)
        {
        	return true;
        }
        else
        {
            printf("I am stuck.\n"); getchar();
            return true;
        }

    }


    bool CavitySignalGenerator::DoGenerate( Signal* aSignal )
    {

        FILE *fp;
        if (fTextFileWriting==1) fp = fopen("incidentfields.txt", "w");


        //n samples for event spacing in Kass.
        int PreEventCounter = 0;

        InitializeFieldPoints(allRxChannels);

        if (!fTransmitter->IsKassiopeia())
        {
        	for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        	{
//        		DriveAntenna(fp, PreEventCounter, index, aSignal, nfilterbins, dtfilter);
        	}  // for loop
        	return true;
        }


        if (fTransmitter->IsKassiopeia())
        {
        	bool fTruth = false;
        	int startingIndex;
            fInterface->fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        	std::thread tKassiopeia (&CavitySignalGenerator::KassiopeiaInit, this, gxml_filename); // spawn new thread

            for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
            {
                if ((!fInterface->fEventInProgress) && (!fInterface->fPreEventInProgress))
                {
                	if (ReceivedKassReady()) fInterface->fPreEventInProgress = true;
                	else
                	{
                		printf("breaking\n");
                		break;
                	}

                	LPROG( lmclog, "LMC ReceivedKassReady" );

                }

                if (fInterface->fPreEventInProgress)  // Locust keeps sampling until Kass event.
                {
                    PreEventCounter += 1;

                    if (((!fTruth)&&(PreEventCounter > fNPreEventSamples))||((fTruth)&&(PreEventCounter > fNPreEventSamples)&&(index%(8192*aSignal->DecimationFactor())==0)  ))// finished pre-samples.  Start event.
                    {
                        fInterface->fPreEventInProgress = false;  // reset.
                        fInterface->fEventInProgress = true;
                        startingIndex = index;
                        LPROG( lmclog, "LMC about to WakeBeforeEvent()" );
                        WakeBeforeEvent();  // trigger Kass event.
                    }
                }

                if (fInterface->fEventInProgress)  // fEventInProgress
                {
                    std::unique_lock< std::mutex >tLock( fInterface->fMutexDigitizer, std::defer_lock );
                	if (!fInterface->fKassEventReady)  // Kass confirms event is underway.
                	{
                        tLock.lock();
                        fInterface->fDigitizerCondition.wait( tLock );
                        if (fInterface->fEventInProgress)
                        {
/*
                      		if (DriveAntenna(fp, startingIndex, index, aSignal, nfilterbins, dtfilter))
                    		{
                                PreEventCounter = 0; // reset
                    		}

                    		else
                    		{
                    			LERROR(lmclog,"The antenna did not respond correctly.  Exiting.\n");
                    			fSkippedSamples = true;
                                tLock.unlock();
                    			break;
                    		}
*/
                        }
                        tLock.unlock();
                	}
                 	else  // diagnose Kass
                 	{
                         tLock.lock();
                         std::this_thread::sleep_for(std::chrono::milliseconds(fThreadCheckTime));
                         if (!fInterface->fKassEventReady)  // Kass event did start.  Continue but skip this sample.
                         {
                         	tLock.unlock();
                         }
                         else  // Kass event has not started.
                         {
                          	if ( fInterface->fEventInProgress )
                          	{
                          		if ( index < fNPreEventSamples+1 ) // Kass never started at all.
                          		{
                         			LERROR(lmclog,"Kass thread is unresponsive.  Exiting.\n");
                             		fKassNeverStarted = true;
                          		}
                             	tLock.unlock(); // Kass either started or not, but is now finished.
                             	break;
                          	}
                          	else  // Kass started an event and quickly terminated it.
                          	{
                         		LWARN(lmclog, "Kass event terminated quickly.\n");
                         		tLock.unlock();
                          	}
                         }
                 	} // diagnose Kass

                } // if fEventInProgress
            }  // for loop

            fInterface->fDoneWithSignalGeneration = true;
            if (fTextFileWriting==1) fclose(fp);
            LPROG( lmclog, "Finished signal loop." );
			fInterface->fWaitBeforeEvent = false;
            WakeBeforeEvent();
            tKassiopeia.join();  // finish thread

            if (fKassNeverStarted == true)
            {
            	throw 2;
            	return false;
            }
            if (fSkippedSamples == true)
            {
            	throw 3;
            	return false;
            }



        }  // fTransmitter->IsKassiopeia()

        return true;

    }

} /* namespace locust */


/*
 * LMCEventHold.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "logger.hh"
#include "LMCEventHold.hh"
#include <csignal>

namespace locust
{

	LOGGER( lmclog, "EventHold" );

    EventHold::EventHold() :
            fTruthOutputFilename("LocustEventProperties.root"),
            fAccumulateTruthInfo( false ),
            fConfigurationComplete( false ),
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    EventHold::EventHold( const EventHold& aOrig ) : KSComponent(),
            fTruthOutputFilename("LocustEventProperties.root"),
            fAccumulateTruthInfo( false ),
            fConfigurationComplete( false ),
            fInterface( aOrig.fInterface )
    {
    }

    EventHold::~EventHold()
    {
    }

    EventHold* EventHold::Clone() const
    {
        return new EventHold( *this );
    }

    bool EventHold::ConfigureByInterface()
    {

        if (!fConfigurationComplete)
        {
            if (fInterface->fConfigureKass)
            {
                const scarab::param_node* aParam = fInterface->fConfigureKass->GetParameters();
                if (!this->Configure( *aParam ))
                {
                    LERROR(lmclog,"Error configuring EventHold class");
                    return false;
                }
            }
            else
            {
                LPROG(lmclog,"EventHold class did not need to be configured.");
                return true;
            }

            OpenFile();

            fConfigurationComplete = true;

        }
        return true;
    }

    bool EventHold::Configure( const scarab::param_node& aParam )
    {
	    if ( aParam.has( "random-track-seed" ) )
	    {
	    	fInterface->anEvent->fRandomSeed = aParam["random-track-seed"]().as_int();
	    }
	    if ( aParam.has( "truth-output-filename" ) )
	    {
	    	fTruthOutputFilename = aParam["truth-output-filename"]().as_string();
	    }
	    if ( aParam.has( "accumulate-truth-info" ) )
	    {
	    	fAccumulateTruthInfo = aParam["accumulate-truth-info"]().as_bool();
	    }


    	return true;
    }



    bool EventHold::OpenEvent()
    {
#ifdef ROOT_FOUND
        fInterface->anEvent = new Event();
        fInterface->anEvent->Initialize();
        fInterface->aTrack = new Track();
        fInterface->aTrack->Initialize();
#endif

        return true;
    }

    bool EventHold::OpenFile()
    {
        std::string tOutputPath = TOSTRING(PB_OUTPUT_DIR);
        std::string sFileName = tOutputPath+"/"+fTruthOutputFilename;
#ifdef ROOT_FOUND
        FileWriter* aRootTreeWriter = RootTreeWriter::get_instance();
        aRootTreeWriter->SetFilename(sFileName);
        if (fAccumulateTruthInfo)
        {
        	// TO-DO:  This option should be used when running pileup.  We will need to
        	// figure out how to explicitly increment the event structure ID, given that the
        	// same (identical) simulation is being run multiple times in this case.
        	aRootTreeWriter->OpenFile("UPDATE");
        }
        else
        {
        	aRootTreeWriter->OpenFile("RECREATE");
        }
        aRootTreeWriter->CloseFile();
#endif
        return true;

    }


    bool EventHold::WriteEvent()
    {
        std::string tOutputPath = TOSTRING(PB_OUTPUT_DIR);
        std::string sFileName = tOutputPath+"/"+fTruthOutputFilename;
#ifdef ROOT_FOUND
        FileWriter* aRootTreeWriter = RootTreeWriter::get_instance();
        aRootTreeWriter->SetFilename(sFileName);
       	aRootTreeWriter->OpenFile("UPDATE");

       	if (fInterface->anEvent->fNTracks == 0) fInterface->anEvent->AddTrack( fInterface->aTrack );
        aRootTreeWriter->WriteEvent( fInterface->anEvent );
        aRootTreeWriter->WriteRunParameters(fInterface->aRunParameter);
        aRootTreeWriter->CloseFile();
#endif
        return true;
    }

    bool EventHold::ExecutePreEventModification(Kassiopeia::KSEvent &anEvent)
    {
    	if ( !ConfigureByInterface() )
    	{
    	    return false;
    	}

        OpenEvent(); // for recording event properties to file.

        LPROG( lmclog, "Kass is waiting for event trigger" );

        fInterface->fDigitizerCondition.notify_one();  // unlock if still locked.
        if(( fInterface->fWaitBeforeEvent ) && (!fInterface->fDoneWithSignalGeneration))
        {
            fInterface->fKassReadyCondition.notify_one();
            std::unique_lock< std::mutex >tLock( fInterface->fMutex );
            fInterface->fPreEventCondition.wait( tLock );
            fInterface->fKassEventReady = false;
            fInterface->fTOld = 0.;  // reset time on event clock
            LPROG( lmclog, "Kass got the event trigger" );
        }
        else
        {
        	printf("Raising sigint to cancel Kassiopeia");
        	raise(SIGINT);
        	return true;
        }
        return true;

    }

    bool EventHold::ExecutePostEventModification(Kassiopeia::KSEvent &anEvent)
    {
        WriteEvent();
        fInterface->fEventInProgress = false;
        fInterface->fDigitizerCondition.notify_one();  // unlock
        LPROG( lmclog, "Kass is waking after event" );
        delete fInterface->anEvent;
        delete fInterface->aTrack;
        return true;
    }

} /* namespace locust */


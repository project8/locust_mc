/*
 * LMCEventHold.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "logger.hh"
#include "LMCEventHold.hh"
#include <csignal>
#include <fstream>


namespace locust
{

	LOGGER( lmclog, "EventHold" );

    EventHold::EventHold() :
            fTruthOutputFilename("LocustEventProperties.root"),
            fAccumulateTruthInfo( false ),
            fConfigurationComplete( false ),
            fEventSeed( 0 ),
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    EventHold::EventHold( const EventHold& aOrig ) : KSComponent(),
            fTruthOutputFilename("LocustEventProperties.root"),
            fAccumulateTruthInfo( false ),
            fConfigurationComplete( false ),
            fEventSeed( 0 ),
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
	    	fEventSeed = aParam["random-track-seed"]().as_int();
	    }
	    else
	    {
	        fEventSeed = -99;
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
        fInterface->anEvent->Initialize( fEventSeed );
        fInterface->aTrack = new Track();
        fInterface->aTrack->Initialize();
#endif

        return true;
    }

    bool EventHold::OpenFile()
    {
        std::string tOutputPath = TOSTRING(PB_OUTPUT_DIR);
        std::string sFileName = tOutputPath+"/"+fTruthOutputFilename;
        fJsonFileName = sFileName;
        const std::string ext(".root");
        if ( sFileName != ext &&
             sFileName.size() > ext.size() &&
             sFileName.substr(sFileName.size() - ext.size()) == ".root" )
        {
            // Remove .root and replace it with .json:
            fJsonFileName = sFileName.substr(0, sFileName.size() - ext.size()) + ".json";
        }
        else
        {
            LERROR(lmclog,"The output file " << fTruthOutputFilename <<"doesn't end in .root");
            exit(-1);
        }
#ifdef ROOT_FOUND
        FileWriter* aRootTreeWriter = RootTreeWriter::get_instance();
        aRootTreeWriter->SetFilename(sFileName);
        if (fAccumulateTruthInfo)
        {
        	aRootTreeWriter->OpenFile("UPDATE");
        }
        else
        {
        	aRootTreeWriter->OpenFile("RECREATE");
        }
        aRootTreeWriter->CloseFile();
#endif

        // Open the json file:
        if (!fAccumulateTruthInfo)
        {
            std::ofstream ost {fJsonFileName, std::ios_base::out};
        }

        return true;

    }


    bool EventHold::WriteJsonFile()
    {
        std::ifstream jsonFile(fJsonFileName); // Open json file for inspection
        bool bNewRun = true;
        std::vector<std::string> v;
        if (jsonFile.is_open())
        {
            std::string line;
            while (std::getline(jsonFile, line))
            {
                LPROG( lmclog, line );
                bNewRun = !line.find("run-id");
                if (line != "}")  // Avoid saving the last "}".  It will be appended below.
                {
                    v.push_back(line);
                }
            }
   	        jsonFile.close();
        }
        else
        {
            LPROG( lmclog, "A json file for meta-data was not found.  Opening a new one now." );
        }

        std::ofstream ost {fJsonFileName, std::ios_base::out};


#ifdef ROOT_FOUND
        if (bNewRun)  // If there are no run parameters in the json file yet, write them now:
        {
            ost << "{\n";
            ost << "    \"run-id\": "<< "\"" << fInterface->aRunParameter->fRunID << "\",\n";
            ost << "    \"run-parameters\": {\n";
            ost << "        \"run-type\": "<< "\"" << fInterface->aRunParameter->fDataType << "\",\n";
            ost << "        \"simulation-type\": "<< "\"" << fInterface->aRunParameter->fSimulationType << "\",\n";
            ost << "        \"simulation-subtype\": "<< "\"" << fInterface->aRunParameter->fSimulationSubType << "\",\n";
            ost << "        \"sampling-freq-mega-hz\": "<< "\"" << fInterface->aRunParameter->fSamplingRateMHz << "\"\n";
            ost << "    },\n";
        }
        else // otherwise re-write the file:
        {
            for (int i = 0; i < v.size(); i++)
            {
                if (i < v.size()-1)
                {
                    ost << v[i] << "\n";
                }
                else
                {
                    ost << v[i] << ",\n";
                }
            }
        }


        // Write the latest event information here:

        ost << "    \"" << fInterface->anEvent->fEventID << "\": {\n";
        ost << "        \"ntracks\": "<< "\"" << fInterface->anEvent->fNTracks << "\",\n";
        for (int i=0; i<fInterface->anEvent->fNTracks; i++)
        {
            ost << "        \"" << fInterface->anEvent->fTrackIDs[i] << "\":\n";
            ost << "         {\n";
            ost << "             \"start-time\": "<< "\"" << fInterface->anEvent->fStartTimes[i] << "\",\n";
            ost << "             \"end-time\": "<< "\"" << fInterface->anEvent->fEndTimes[i] << "\",\n";
            ost << "             \"output-avg-frequency\": "<< "\"" << fInterface->anEvent->fOutputAvgFrequencies[i] << "\",\n";
            ost << "             \"pitch-angle\": "<< "\"" << fInterface->anEvent->fPitchAngles[i] << "\",\n";
            ost << "             \"avg-axial-frequency\": "<< "\"" << fInterface->anEvent->fAvgAxialFrequencies[i] << "\"\n";
            if (i < fInterface->anEvent->fNTracks-1)
            {
                ost << "         },\n";
            }
            else
            {
            	ost << "         }\n";
            }
        }
        ost << "    }\n";
        ost << "}\n";

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
        WriteJsonFile();
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
#ifdef ROOT_FOUND
        delete fInterface->anEvent;
        delete fInterface->aTrack;
        delete fInterface->aRunParameter;
#endif
        return true;
    }

} /* namespace locust */


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
			fMetadataTag("none"),
            fAccumulateTruthInfo( false ),
            fConfigurationComplete( false ),
            fEventSeed( 0 ),
            fConfiguredEMin( 0. ),
            fConfiguredPitchMin( 0. ),
            fConfiguredXMin( 0. ),
            fConfiguredYMin( 0. ),
            fConfiguredZMin( 0. ),
            fEventCounter ( 0 ),
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    EventHold::EventHold( const EventHold& aOrig ) : KSComponent(),
            fTruthOutputFilename("LocustEventProperties.root"),
            fAccumulateTruthInfo( false ),
            fConfigurationComplete( false ),
            fEventSeed( 0 ),
            fConfiguredEMin( 0. ),
            fConfiguredPitchMin( 0. ),
            fConfiguredXMin( 0. ),
            fConfiguredYMin( 0. ),
            fConfiguredZMin( 0. ),
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
        if ( aParam.has( "metadata-tag" ) )
        {
            fMetadataTag = aParam["metadata-tag"]().as_string();
        }
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
	    if ( aParam.has( "ks-starting-energy-min" ) )
	    {
            fConfiguredEMin = aParam["ks-starting-energy-min"]().as_double();
	    }
	    if ( aParam.has( "ks-starting-pitch-min" ) )
	    {
            fConfiguredPitchMin = aParam["ks-starting-pitch-min"]().as_double();
	    }
	    if ( aParam.has( "ks-starting-xpos-min" ) )
	    {
            fConfiguredXMin = aParam["ks-starting-xpos-min"]().as_double();
	    }
	    if ( aParam.has( "ks-starting-ypos-min" ) )
	    {
            fConfiguredYMin = aParam["ks-starting-ypos-min"]().as_double();
	    }
	    if ( aParam.has( "ks-starting-zpos-min" ) )
	    {
            fConfiguredZMin = aParam["ks-starting-zpos-min"]().as_double();
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
        fEventCounter = 0;
        std::vector<std::string> v;
        if (jsonFile.is_open())
        {
            std::string line;
            while (std::getline(jsonFile, line))
            {
                bNewRun = !line.find("run-id");
                //increment the event counter
                if ( line.find("\"event-tag\"") != std::string::npos ) fEventCounter += 1;
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

        FILE *file = std::fopen(fJsonFileName.c_str(), "w");


#ifdef ROOT_FOUND
        if (bNewRun)  // If there are no run parameters in the json file yet, write them now:
        {
            struct timeval tv;
            gettimeofday(&tv, NULL);
            int tMillisec = int( tv.tv_usec / 1000);

            fprintf(file, "{\n");
            fprintf(file, "    \"run-id\": \"%ld-%03d\",\n", fInterface->aRunParameter->fRunID, tMillisec);
            fprintf(file, "    \"run-parameters\": {\n");
            fprintf(file, "        \"run-type\": \"%s\",\n", fInterface->aRunParameter->fDataType.c_str());
            fprintf(file, "        \"simulation-type\": \"%s\",\n", fInterface->aRunParameter->fSimulationType.c_str());
            fprintf(file, "        \"simulation-subtype\": \"%s\",\n", fInterface->aRunParameter->fSimulationSubType.c_str());
            fprintf(file, "        \"user-defined-tag\": \"%s\",\n", fMetadataTag.c_str());
            fprintf(file, "        \"sampling-freq-mega-hz\": \"%.1f\",\n", fInterface->aRunParameter->fSamplingRateMHz);
            fprintf(file, "        \"configured-e-min\": \"%12.6f\",\n", fConfiguredEMin);
            fprintf(file, "        \"configured-pitch-min\": \"%12.9f\",\n", fConfiguredPitchMin);
            fprintf(file, "        \"configured-x-min\": \"%11.9f\",\n", fConfiguredXMin);
            fprintf(file, "        \"configured-y-min\": \"%11.9f\",\n", fConfiguredYMin);
            fprintf(file, "        \"configured-z-min\": \"%11.9f\"\n", fConfiguredZMin);
            fprintf(file, "    },\n");
            fprintf(file, "    \"nevents\": %d,\n", fEventCounter+1);
        }
        else // otherwise re-write the file:
        {
            for (int i = 0; i < v.size(); i++)
            {
                if (i < v.size()-1)
                {
                    if ( v[i].find("\"nevents\"") != std::string::npos )
                    {
                    	fprintf(file, "    \"nevents\": %d,\n", fEventCounter+1);
                    }
                    else
                    {
                        fprintf(file,"%s\n", v[i].c_str());
                    }
                }
                else
                {
                    fprintf(file,"%s,\n", v[i].c_str());
                }
            }
        }


        // Write the latest event information here:

        fprintf(file,"    \"%d\": {\n", fEventCounter);
        fprintf(file,"        \"event-tag\": \"%ld\",\n", fInterface->anEvent->fEventID);
        fprintf(file,"        \"kassiopeia-seed\": %ld,\n", fInterface->aRunParameter->fKassiopeiaSeed);
        fprintf(file,"        \"track-length-seed\": %d,\n", fInterface->aRunParameter->fTrackLengthSeed);
        fprintf(file,"        \"track-delay-seed\": %d,\n", fInterface->aRunParameter->fTrackDelaySeed);
        fprintf(file,"        \"ntracks\": %d,\n", fInterface->anEvent->fNTracks);
        for (int i=0; i<fInterface->anEvent->fNTracks; i++)
        {
            fprintf(file,"        \"%d\":\n", fInterface->anEvent->fTrackIDs[i]);
            fprintf(file,"         {\n");
            fprintf(file,"             \"start-time\": %g,\n", fInterface->anEvent->fStartTimes[i]);
            fprintf(file,"             \"end-time\": %g,\n", fInterface->anEvent->fEndTimes[i]);
            fprintf(file,"             \"energy-ev\": %12.10g,\n", fInterface->anEvent->fStartingEnergies_eV[i]);
            fprintf(file,"             \"start-radius\": %12.10g,\n", fInterface->anEvent->fRadii[i]);
            fprintf(file,"             \"start-radial-phase\": %12.10g,\n", fInterface->anEvent->fRadialPhases[i]);
            fprintf(file,"             \"start-guiding-center-x\": %12.10g,\n", fInterface->anEvent->fStartGuidingCentersX[i]);
            fprintf(file,"             \"start-guiding-center-y\": %12.10g,\n", fInterface->anEvent->fStartGuidingCentersY[i]);
            fprintf(file,"             \"start-guiding-center-z\": %12.10g,\n", fInterface->anEvent->fStartGuidingCentersZ[i]);
            fprintf(file,"             \"output-avg-frequency\": %12.10g,\n", fInterface->anEvent->fOutputAvgFrequencies[i]);
            fprintf(file,"             \"output-track-start-frequency\": %12.10g,\n", fInterface->anEvent->fTrackOutputStartFrequencies[i]);
            fprintf(file,"             \"output-inst-start-frequency\": %12.10g,\n", fInterface->anEvent->fOutputStartFrequencies[i]);
            fprintf(file,"             \"pitch-angle\": %.6f,\n", fInterface->anEvent->fPitchAngles[i]);
            fprintf(file,"             \"slope\": %g,\n", fInterface->anEvent->fSlopes[i]);
            fprintf(file,"             \"avg-axial-frequency\": %g\n", fInterface->anEvent->fAvgAxialFrequencies[i]);
            if (i < fInterface->anEvent->fNTracks-1)
            {
                fprintf(file,"         },\n");
            }
            else
            {
            	fprintf(file,"         }\n");
            }
        }
        fprintf(file,"    }\n");
        fprintf(file,"}\n");

#endif
        std::fclose(file);

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
#endif
        return true;
    }

} /* namespace locust */


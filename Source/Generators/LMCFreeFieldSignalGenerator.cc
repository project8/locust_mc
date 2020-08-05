/*
 * LMCFreeFieldSignalGenerator.cc
 *
 *  Created on: Mar 12, 2017
 *      Author: buzinsky 
 */

#include "LMCFreeFieldSignalGenerator.hh"
#include "LMCRunKassiopeia.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

//#include "LMCGlobalsDeclaration.hh"
#include "LMCHFSSReader.hh"

namespace locust
{
    LOGGER( lmclog, "FreeFieldSignalGenerator" );

    MT_REGISTER_GENERATOR(FreeFieldSignalGenerator, "freefield-signal");

    FreeFieldSignalGenerator::FreeFieldSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fLO_Frequency( 0.),
		fNPreEventSamples( 150000 ),
        fArrayRadius(0.),
        fElementSpacing(0.),
        fNElementsPerStrip(0.),
        fCorporateFeed(1),
        fPileupSeed( 0.),
        fPileupMode( false),
        gxml_filename("blank.xml"),
		fInterface( new KassLocustInterface() )

    {

    	fRequiredSignalState = Signal::kTime;

        KLInterfaceBootstrapper::get_instance()->SetInterface( fInterface );

    }

    FreeFieldSignalGenerator::~FreeFieldSignalGenerator()
    {
    }

    bool FreeFieldSignalGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "xml-filename" ) )
        {
            gxml_filename = aParam["xml-filename"]().as_string();
        }

        if( aParam.has( "array-radius" ) )
        {
            fArrayRadius = aParam["array-radius"]().as_double();
        }

        if( aParam.has( "nelements-per-strip" ) )
        {
            fNElementsPerStrip = aParam["nelements-per-strip"]().as_double();
        }

        if( aParam.has( "element-spacing" ) )
        {
            fElementSpacing = aParam["element-spacing"]().as_double();
        }

        if( aParam.has( "pileup" ) )
        {
            fPileupMode = aParam["pileup"]().as_bool();
        }

        if( aParam.has( "pileup-seed" ) )
        {
            fPileupSeed = aParam["pileup-seed"]().as_int();
        }
        if( aParam.has( "feed" ) )
        {
            std::string feedInput = aParam["feed"]().as_string();
            if(feedInput == "series")
                fCorporateFeed = false;
        }


        return true;
    }

    void FreeFieldSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    void FreeFieldSignalGenerator::KassiopeiaInit(const std::string &aFile)
    {
        RunKassiopeia tRunKassiopeia;
        tRunKassiopeia.Run(aFile, fInterface);
        return;
    }



    bool FreeFieldSignalGenerator::WakeBeforeEvent()
    {
        fInterface->fPreEventCondition.notify_one();
        return true;
    }

    bool FreeFieldSignalGenerator::ReceivedKassReady()
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        printf("LMC about to wait ..\n ");

        if( ! fInterface->fKassEventReady)
        {
            std::unique_lock< std::mutex >tLock( fInterface->fKassReadyMutex );
            fInterface->fKassReadyCondition.wait( tLock );
        }

        if (fInterface->fFalseStartKassiopeia)  // workaround for some Macs
        {
            std::unique_lock< std::mutex >tLock( fInterface->fKassReadyMutex );
            fInterface->fKassReadyCondition.wait( tLock );

        }

        return true;
    }



    void* FreeFieldSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal)
    {

        const int signalSize = aSignal->TimeSize();

        //Receiver Properties
        double tReceiverTime = fInterface->fTOld;
        double tRetardedTime = 0.; //Retarded time of particle corresponding to when emission occurs, reaching receiver at tReceiverTime
        const double tLocustStep = 1. / (fAcquisitionRate * 1e6 * aSignal->DecimationFactor());

        PatchAntenna *currentPatch;
        unsigned tTotalPatchIndex = 0;

        fFieldSolver.SetParticleHistory(fInterface->fParticleHistory);

        for(int channelIndex = 0; channelIndex < allChannels.size(); ++channelIndex)
        {
            for(int patchIndex = 0; patchIndex < allChannels[channelIndex].size(); ++patchIndex)
            {
                currentPatch = &allChannels[channelIndex][patchIndex]; 

                fFieldSolver.SetFieldEvent(tReceiverTime, tTotalPatchIndex);
                fFieldSolver.SolveFieldSolutions();

                LMCThreeVector tRadiatedElectricField = fFieldSolver.GetElectricField();
                LMCThreeVector tRadiatedMagneticField = fFieldSolver.GetMagneticField();
                locust::Particle tCurrentParticle = fFieldSolver.GetRetardedParticle();

                //////////////////////////////////////////////
                LMCThreeVector tDirection = currentPatch->GetPosition() - tCurrentParticle.GetPosition(true);
                double tVelZ = tCurrentParticle.GetVelocity(true).Z();
                double tCosTheta =  tVelZ * tDirection.Z() /  tDirection.Magnitude() / fabs(tVelZ);
                double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency() / ( 1. - fabs(tVelZ) / LMCConst::C() * tCosTheta);

                currentPatch->SetInstantaneousFrequency( tDopplerFrequency / (2. * LMCConst::Pi() ));
                currentPatch->SetIncidentElectricField( tRadiatedElectricField );
                currentPatch->SetIncidentMagneticField( tRadiatedMagneticField );

                aSignal->LongSignalTimeComplex()[channelIndex*signalSize*aSignal->DecimationFactor() + index][0] += currentPatch->GetVoltage();
                //aSignal->LongSignalTimeComplex()[channelIndex*signalSize*aSignal->DecimationFactor() + index][1] += currentPatch->GetVoltage();

                ++tTotalPatchIndex;

            } // z_position waveguide element stepping loop.
        } // nChannels loop.

        fInterface->fTOld += tLocustStep;

        return 0;
    }


    void FreeFieldSignalGenerator::InitializePatchArray()
    {

        const unsigned nChannels = fNChannels;
        const int nReceivers = fNElementsPerStrip; //Number of receivers per channel

        const double patchSpacingZ = fElementSpacing;
        const double patchRadius = fArrayRadius;
        double zPosition;
        double theta;
        const double dThetaArray = 2. * LMCConst::Pi() / nChannels; //Divide the circle into nChannels

        PatchAntenna modelPatch;

        allChannels.resize(nChannels);

        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            theta = channelIndex * dThetaArray;

            for(int receiverIndex = 0; receiverIndex < nReceivers; ++receiverIndex)
            {
                zPosition =  (receiverIndex - (nReceivers - 1.) /2.) * patchSpacingZ;

                modelPatch.SetCenterPosition({patchRadius * cos(theta) , patchRadius * sin(theta) , zPosition }); 
                modelPatch.SetPolarizationDirection({sin(theta), -cos(theta), 0.}); 
                modelPatch.SetNormalDirection({-cos(theta), -sin(theta), 0.}); //Say normals point inwards
                allChannels[channelIndex].AddReceiver(modelPatch);
                fFieldSolver.AddFieldPoint(modelPatch.GetPosition());
            }
        }
    }

    bool FreeFieldSignalGenerator::DoGenerate( Signal* aSignal )
    {
        InitializePatchArray();

        // Initialize random number generator for pileup
        if(fPileupMode)
        { 
            if(fPileupSeed!=0)
                srand(fPileupSeed);
            else 
                srand(time(NULL));
        }

        //n samples for event spacing.
        int PreEventCounter = 0;
        fInterface->fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        std::thread tKassiopeia (&FreeFieldSignalGenerator::KassiopeiaInit, this, gxml_filename);     // spawn new thread

        unsigned index = 0;

        while( index < aSignal->DecimationFactor()*aSignal->TimeSize())
        {
            if ((!fInterface->fEventInProgress) && (!fInterface->fPreEventInProgress))
            {
            	if (ReceivedKassReady()) fInterface->fPreEventInProgress = true;
            	fInterface->fPreEventInProgress = true;
            	printf("LMC says it ReceivedKassReady()\n");

            }

            if (fInterface->fPreEventInProgress)
            {
                PreEventCounter += 1;

                if (PreEventCounter > fNPreEventSamples) // finished pre-samples.  Start event.
                {
                    fInterface->fPreEventInProgress = false;  // reset.
                    fInterface->fEventInProgress = true;
                    printf("LMC about to WakeBeforeEvent()\n");
                    WakeBeforeEvent();  // trigger Kass event.
                    if(fPileupMode)
                      {
                          index = rand() % aSignal->DecimationFactor()*aSignal->TimeSize();
                      }
                }
            }

            if (fInterface->fEventInProgress)  // fEventInProgress
            {
                    std::unique_lock< std::mutex >tLock( fInterface->fMutexDigitizer, std::defer_lock );
                    tLock.lock();
                    fInterface->fDigitizerCondition.wait( tLock );
                    if (fInterface->fEventInProgress)
                    {
                        DriveAntenna(PreEventCounter, index, aSignal);
                        PreEventCounter = 0; // reset
                    }
                    tLock.unlock();
            }
        }  // for loop


        fInterface->fDoneWithSignalGeneration = true;  // tell LMCCyclotronRadExtractor
        WakeBeforeEvent();  // trigger one last Kass event if we are locked up.
        tKassiopeia.join();  // finish thread



        return true;
    }

} /* namespace locust */

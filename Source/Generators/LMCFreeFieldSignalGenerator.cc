/*
 * LMCFreeFieldSignalGenerator.cc
 *
 *  Created on: Mar 12, 2017
 *      Author: buzinsky 
 */

#include "LMCFreeFieldSignalGenerator.hh"
#include "LMCEventHold.hh"
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
        //fWriteNFD(0.),
        fLO_Frequency( 0.),
        fArrayRadius(0.),
        fPatchSpacing(0.),
        fNPatchesPerStrip(0.),
        fCorporateFeed(1),
        fPileupSeed( 0.),
        fPileupMode( false),
        gxml_filename("blank.xml")
    {
        fRequiredSignalState = Signal::kTime;
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

        if( aParam.has( "npatches-per-strip" ) )
        {
            fNPatchesPerStrip = aParam["npatches-per-strip"]().as_double();
        }

        if( aParam.has( "patch-spacing" ) )
        {
            fPatchSpacing = aParam["patch-spacing"]().as_double();
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


    static void* KassiopeiaInit(const std::string &aFile)
    {
        RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
        RunKassiopeia1->Run(aFile);
        delete RunKassiopeia1;

        return 0;
    }

    static void WakeBeforeEvent()
    {
        fPreEventCondition.notify_one();
        return;
    }

    static bool ReceivedKassReady()
    {
        printf("LMC about to wait ..\n");

        std::unique_lock< std::mutex >tLock( fKassReadyMutex);
        fKassReadyCondition.wait( tLock, [](){return fKassEventReady;} );
        printf("LMC Got the fKassReadyCondition signal\n");

        return true;
    }

    double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition )
    {
        return aReceiverTime - aParticleTime - (aReceiverPosition - aParticlePosition).Magnitude() / LMCConst::C();
    }

    double GetStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval)
    {
        double tRetardedTime = aParticle.GetTime(true);
        return tRetardedTime + aSpaceTimeInterval;
    }


    void* FreeFieldSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal)
    {

        locust::Particle tCurrentParticle = fParticleHistory.back();
        int CurrentIndex;
        const int signalSize = aSignal->TimeSize();

        const double kassiopeiaTimeStep = fabs(fParticleHistory[0].GetTime() - fParticleHistory[1].GetTime());
        const int historySize = fParticleHistory.size();

        //Receiver Properties
        double tReceiverTime = t_old;

        double tRetardedTime = 0.; //Retarded time of particle corresponding to when emission occurs, reaching receiver at tReceiverTime

        double tSpaceTimeInterval=99.;
        double dtRetarded=0;
        double tTolerance=1e-23;

        PatchAntenna *currentPatch;

        for(int channelIndex = 0; channelIndex < allChannels.size(); ++channelIndex)
        {
            for(int patchIndex = 0; patchIndex < allChannels[channelIndex].size(); ++patchIndex)
            {
                currentPatch = &allChannels[channelIndex][patchIndex]; 
                //tReceiverTime = t_old - fabs(currentPatch->GetPosition().Z()) / LMCConst::C();

                if(fParticleHistory.front().GetTime()<=3.*kassiopeiaTimeStep)
                {
                    fParticleHistory.front().Interpolate(0);
                    if(GetSpaceTimeInterval(fParticleHistory.front().GetTime(true), tReceiverTime , fParticleHistory.front().GetPosition(true), currentPatch->GetPosition() ) < 0 )
                    {
                        //printf("Skipping! out of Bounds!: tReceiverTime=%e\n",tReceiverTime);
                        continue;
                    }
                }

                if(currentPatch->GetPreviousRetardedIndex() == -99.)
                {
                    CurrentIndex=FindNode(tReceiverTime);
                    tCurrentParticle = fParticleHistory[CurrentIndex];
                    const double vGroup = 1.674e8;
                    tRetardedTime = tReceiverTime - (tCurrentParticle.GetPosition() - currentPatch->GetPosition() ).Magnitude() /  vGroup;
                    if(tRetardedTime < 0) 
                    {
                        tRetardedTime = 0;
                        CurrentIndex = 0;
                    }
                }
                else
                {
                    CurrentIndex = currentPatch->GetPreviousRetardedIndex();
                    tRetardedTime = currentPatch->GetPreviousRetardedTime() + 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
                }


                CurrentIndex = FindNode(tRetardedTime);

                tCurrentParticle = fParticleHistory[CurrentIndex];
                tCurrentParticle.Interpolate(tRetardedTime);
                tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), currentPatch->GetPosition());

                double tOldSpaceTimeInterval=99.;

                //Converge to root
                for(int j=0;j<25;++j)
                {
                    //++tAverageIterations;

                    tRetardedTime = GetStepRoot(tCurrentParticle, tReceiverTime, currentPatch->GetPosition(), tSpaceTimeInterval);
                    tCurrentParticle.Interpolate(tRetardedTime);

                    //Change the kassiopeia step we expand around if the interpolation time displacement is too large
                    if(fabs(tCurrentParticle.GetTime(true) - tCurrentParticle.GetTime(false)) > kassiopeiaTimeStep)
                    {
                        CurrentIndex=FindNode(tRetardedTime);
                        tCurrentParticle=fParticleHistory[CurrentIndex];
                        tCurrentParticle.Interpolate(tRetardedTime);
                    }

                    tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), currentPatch->GetPosition());
                    tOldSpaceTimeInterval = tSpaceTimeInterval;
                }

                currentPatch->SetPreviousRetardedIndex(CurrentIndex);
                currentPatch->SetPreviousRetardedTime(tRetardedTime);

                LMCThreeVector tDirection = currentPatch->GetPosition() - tCurrentParticle.GetPosition(true);
                double tVelZ = tCurrentParticle.GetVelocity(true).Z();
                double tCosTheta =  tVelZ * tDirection.Z() /  tDirection.Magnitude() / fabs(tVelZ);
                double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency() / ( 1. - fabs(tVelZ) / LMCConst::C() * tCosTheta);

                currentPatch->SetInstantaneousFrequency( tDopplerFrequency / (2. * LMCConst::Pi() ));
                currentPatch->SetIncidentElectricField( tCurrentParticle.CalculateElectricField(currentPatch->GetPosition() ));
                currentPatch->SetIncidentMagneticField( tCurrentParticle.CalculateMagneticField(currentPatch->GetPosition() ));

                aSignal->LongSignalTimeComplex()[channelIndex*signalSize*aSignal->DecimationFactor() + index][0] += currentPatch->GetVoltage();
                //aSignal->LongSignalTimeComplex()[channelIndex*signalSize*aSignal->DecimationFactor() + index][1] += currentPatch->GetVoltage();



            } // z_position waveguide element stepping loop.
        } // nChannels loop.

        t_old += 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        return 0;
    }


    //Return index of fParticleHistory particle closest to the time we are evaluating
    int FreeFieldSignalGenerator::FindNode(double tNew) const
    {
        std::deque<locust::Particle>::iterator it;

        //Get iterator pointing to particle step closest to tNew
        it = std::upper_bound( fParticleHistory.begin() , fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

        int tNodeIndex = it - fParticleHistory.begin();

        return tNodeIndex;
    }

    void FreeFieldSignalGenerator::InitializePatchArray()
    {

        const unsigned nChannels = fNChannels;
        const int nReceivers = fNPatchesPerStrip; //Number of receivers per channel

        const double patchSpacingZ = fPatchSpacing;
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
        const int NPreEventSamples = 150000;
        fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        std::thread Kassiopeia (KassiopeiaInit, gxml_filename);     // spawn new thread

        fRunInProgress = true;
        unsigned index = 0;

        while( index < aSignal->DecimationFactor()*aSignal->TimeSize() && fRunInProgress)
        {
            if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
            {
                if (ReceivedKassReady()) fPreEventInProgress = true;
            }

            if (fPreEventInProgress)
            {
                PreEventCounter += 1;
                if (PreEventCounter > NPreEventSamples)  // finished noise samples.  Start event.
                {
                    fPreEventInProgress = false;  // reset.
                    fEventInProgress = true;
                    WakeBeforeEvent();  // trigger Kass event.
                    if(fPileupMode)
                    {
                        index = rand() % aSignal->DecimationFactor()*aSignal->TimeSize();
                    }

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
                        DriveAntenna(PreEventCounter, index, aSignal);

                        PreEventCounter = 0; // reset
                    }
                    tLock.unlock();
                }

            ++index;

        }  // for loop

        //if(fWriteNFD) NFDWrite();

        //delete [] ImaginarySignal;

        // trigger any remaining events in Kassiopeia so that its thread can finish.
        while (fRunInProgress)
        {
            if (fRunInProgress)
                if (ReceivedKassReady()) WakeBeforeEvent();
        }

        Kassiopeia.join();  // wait for Kassiopeia to finish.

        return true;
    }

} /* namespace locust */

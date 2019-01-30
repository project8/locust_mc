/*
 * LMCPatchSignalGenerator.cc
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
 */

#include "LMCPatchSignalGenerator.hh"
#include "LMCEventHold.hh"
#include "LMCRunKassiopeia.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>

#include "LMCGlobalsDeclaration.hh"
#include "LMCDigitizer.hh"
#include <chrono>


namespace locust
{
    LOGGER( lmclog, "PatchSignalGenerator" );

    MT_REGISTER_GENERATOR(PatchSignalGenerator, "patch-signal");

    PatchSignalGenerator::PatchSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fLO_Frequency( 0.),
        fArrayRadius( 0. ),
        fNPatchesPerStrip( 0. ),
        fPatchSpacing( 0. ),
        fPowerCombiner( 0 ),
        gxml_filename("blank.xml"),
        phiLO_t(0.),
        VoltagePhase_t {0.}
    {
        fRequiredSignalState = Signal::kTime;
    }

    PatchSignalGenerator::~PatchSignalGenerator()
    {
    }

    bool PatchSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        if( aParam->has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam->get_value< double >( "lo-frequency" );
        }
        if( aParam->has( "array-radius" ) )
        {
            fArrayRadius = aParam->get_value< double >( "array-radius" );
        }
        if( aParam->has( "npatches-per-strip" ) )
        {
            fNPatchesPerStrip = aParam->get_value< int >( "npatches-per-strip" );
        }
        if( aParam->has( "patch-spacing" ) )
        {
            fPatchSpacing = aParam->get_value< double >( "patch-spacing" );
        }
        if( aParam->has( "xml-filename" ) )
        {
            gxml_filename = aParam->get_value< std::string >( "xml-filename" );
        }
        if( aParam->has( "feed" ) )
        {
            if (aParam->get_value< std::string >( "feed" ) == "corporate")
                fPowerCombiner = 0;  // default
            else if (aParam->get_value< std::string >( "feed" ) == "series")
                fPowerCombiner = 1;
            else if (aParam->get_value< std::string >( "feed" ) == "quadraturefeed")
                fPowerCombiner = 2;
            else
            	fPowerCombiner = 0;  // default
        }

        return true;
    }

    void PatchSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    static void* KassiopeiaInit(const std::string &aFile)
    {
        //RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
        RunKassiopeia RunKassiopeia1;
        RunKassiopeia1.Run(aFile);
        RunKassiopeia1.~RunKassiopeia();
        //delete RunKassiopeia1;

        return 0;
    }



    static void WakeBeforeEvent()
    {
        fPreEventCondition.notify_one();
        return;
    }

    static bool ReceivedKassReady()
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



    double PatchSignalGenerator::GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition )
    {
        return aReceiverTime - aParticleTime - (aReceiverPosition - aParticlePosition).Magnitude() / LMCConst::C();
    }


    double GetPatchStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval)
    {
        double tRetardedTime = aParticle.GetTime(true);
        return tRetardedTime + aSpaceTimeInterval;
    }


    double GetMismatchFactor(double f)  
    {
        //f /= 2.*LMCConst::Pi();
        // placeholder = 1 - mag(S11)
        // fit to HFSS output
        //double MismatchFactor = 1. - (-5.39e16 / ((f-25.9141e9)*(f-25.9141e9) + 7.23e16) + 0.88);
        //    printf("dopplerfrequency is %f and mismatchfactor is %g\n", f, MismatchFactor);  getchar();
        double MismatchFactor = 0.85;  // punt.
        return MismatchFactor;
    }

    double GetAOIFactor(LMCThreeVector IncidentKVector, double PatchPhi)
    {
        LMCThreeVector PatchNormalVector;
        PatchNormalVector.SetComponents(cos(PatchPhi), sin(PatchPhi), 0.0);
        double AOIFactor = fabs(IncidentKVector.Unit().Dot(PatchNormalVector));
        //printf("cos aoi is %f\n", AOIFactor);
        return AOIFactor;
    }

    double GetVoltageAmpFromPlaneWave()
    {
        double AntennaFactor = 1./420.;

        // S = epsilon0 c E0^2 / 2.  // power/area
        //  0.6e-21 W/Hz * 24.e3 Hz / (0.00375*0.002916) = S = 1.3e-12 W/m^2
        // We should detect 0.6e-21 W/Hz * 24.e3 Hz in Katydid.
        // E0 = sqrt(2.*S/epsilon0/c)
        // effective patch area 0.00004583662 m^2 

        double S = 0.6e-21*24.e3/(0.00004271);  // W/m^2, effective aperture.
        double E0 = sqrt(2.*S/(LMCConst::EpsNull() * LMCConst::C()));
        //    double E0 = 1.0; // V/m, test case
        double amplitude = E0*AntennaFactor;  // volts
        return amplitude;


    }


    // voltage amplitude induced at patch.
    double GetVoltageAmplitude(LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double PatchPhi, double DopplerFrequency)
    {
        double AntennaFactor = 1./400.;
        double MismatchFactor = GetMismatchFactor(DopplerFrequency);
        double AOIFactor = GetAOIFactor(IncidentKVector, PatchPhi);  // k dot patchnormal
        LMCThreeVector PatchPolarizationVector;
        PatchPolarizationVector.SetComponents(-sin(PatchPhi), cos(PatchPhi), 0.0);
        double VoltageAmplitude = fabs( AntennaFactor * IncidentElectricField.Dot(PatchPolarizationVector) * MismatchFactor * AOIFactor);
        //double VoltageAmplitude = fabs( AntennaFactor * IncidentElectricField.Magnitude()); // test case.  

        //    if (VoltageAmplitude>0.) {printf("IncidentElectricField.Dot(PatchPolarizationVector) is %g and VoltageAmplitude is %g\n", IncidentElectricField.Dot(PatchPolarizationVector), VoltageAmplitude); getchar();}
        return VoltageAmplitude;
    }


    // z-index ranges from 0 to npatches-per-strip-1.
    void PatchSignalGenerator::AddOnePatchVoltageToStripSum(Signal* aSignal, double VoltageAmplitude, double VoltagePhase, double phi_LO, unsigned channelindex, unsigned z_index, double DopplerFrequency)
    {
    	PowerCombiner aPowerCombiner;
        if (fPowerCombiner == 1)  // series feed
        {
        	//lossless series feed with amp at one end:
        	VoltagePhase += aPowerCombiner.GetLinePhaseCorr(z_index, DopplerFrequency);
        }
        if (fPowerCombiner == 2) // quadrature feed
        {
        	// assume 2PI delay between junctions, so we don't calculated phase mismatches.
        	// instead calculate damping on voltage amplitude:
        	int njunctions = (int)fabs(z_index - fNPatchesPerStrip/2);
            VoltageAmplitude *= aPowerCombiner.GetVoltageDamping(njunctions);
        }

        //if (VoltageAmplitude>0.) {printf("voltageamplitude is %g\n", VoltageAmplitude); getchar();}
        aSignal->LongSignalTimeComplex()[channelindex][0] += VoltageAmplitude * cos(VoltagePhase - phi_LO);
        aSignal->LongSignalTimeComplex()[channelindex][1] += VoltageAmplitude * sin(VoltagePhase - phi_LO);

    }


    void* PatchSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal)
    {
        if (PreEventCounter > 0)  // new event starting.                                                    
        {
            // initialize patch voltage phases.                                                             
            for (unsigned i=0; i < sizeof(VoltagePhase_t)/sizeof(VoltagePhase_t[0]); i++)
            {
                VoltagePhase_t[i] = {0.};
            }
        }

        locust::Particle tCurrentParticle = fParticleHistory.back();
        int CurrentIndex;
        const int signalSize = aSignal->TimeSize();

        const double kassiopeiaTimeStep = fabs(fParticleHistory[0].GetTime() - fParticleHistory[1].GetTime());
        const int historySize = fParticleHistory.size();
        unsigned sampleIndex = 0;

        //Receiver Properties
        phiLO_t += 2. * LMCConst::Pi() * fLO_Frequency * fDigitizerTimeStep; 
        double tReceiverTime = t_old;
        double tRetardedTime = 0.; //Retarded time of particle corresponding to when emission occurs, reaching receiver at tReceiverTime

        double tSpaceTimeInterval=99.;
        double dtRetarded=0;
        double tTolerance=1e-23;

        PatchAntenna *currentPatch;

        for(int channelIndex = 0; channelIndex < allChannels.size(); ++channelIndex)
        {
            double PatchPhi = (double)channelIndex*360./allChannels.size()*LMCConst::Pi()/180.; // radians.    
            for(int patchIndex = 0; patchIndex < allChannels[channelIndex].size(); ++patchIndex)
            {
                currentPatch = &allChannels[channelIndex][patchIndex]; 
                sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample

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
                    tRetardedTime = tReceiverTime - (tCurrentParticle.GetPosition() - currentPatch->GetPosition() ).Magnitude() /  LMCConst::C();
                }
                else
                {
                    CurrentIndex = currentPatch->GetPreviousRetardedIndex();
                    tRetardedTime = currentPatch->GetPreviousRetardedTime() + fDigitizerTimeStep;
                }

                CurrentIndex = FindNode(tRetardedTime);
                CurrentIndex = std::min(std::max(CurrentIndex,0) , historySize - 1);

                tCurrentParticle = fParticleHistory[CurrentIndex];
                tCurrentParticle.Interpolate(tRetardedTime);
                tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), currentPatch->GetPosition());

                double tOldSpaceTimeInterval=99.;

                //Converge to root
                for(int j=0;j<25;++j)
                {
                    tRetardedTime = GetPatchStepRoot(tCurrentParticle, tReceiverTime, currentPatch->GetPosition(), tSpaceTimeInterval);
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


                if (VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex]>0.)  // not first sample                                                        
                {
                    VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex] += tDopplerFrequency * fDigitizerTimeStep;
                }
                else  // if this is the first light at this patch, the voltage phase doesn't advance for the full dt.
                {              
                    VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex] += tDopplerFrequency * tRetardedTime;
                    //printf("tDopplerFrequency is %g\n", tDopplerFrequency); getchar();
                }      

                double tVoltageAmplitude = GetVoltageAmplitude(tCurrentParticle.CalculateElectricField(currentPatch->GetPosition()), tCurrentParticle.CalculateElectricField(currentPatch->GetPosition()).Cross(tCurrentParticle.CalculateMagneticField(currentPatch->GetPosition())), PatchPhi, tDopplerFrequency);
                //printf("tVoltageAmplitude is %g\n", tVoltageAmplitude); getchar();

                AddOnePatchVoltageToStripSum(aSignal, tVoltageAmplitude, VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex], phiLO_t, sampleIndex, patchIndex, tDopplerFrequency);

            } // patch loop
            //printf("channel %d voltage is %g\n", channelIndex, aSignal->LongSignalTimeComplex()[sampleIndex][0]); getchar();

        } // channels loop

        t_old += fDigitizerTimeStep;

        return 0;
    }

    //Return index of fParticleHistory particle closest to the time we are evaluating
    int PatchSignalGenerator::FindNode(double tNew) const
    {
        std::deque<locust::Particle>::iterator it;
        it = std::upper_bound( fParticleHistory.begin() , fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

        int tNodeIndex = it - fParticleHistory.begin();

        return tNodeIndex;
    }

    void PatchSignalGenerator::InitializePatchArray()
    {

        const unsigned nChannels = fNChannels;
        const int nReceivers = fNPatchesPerStrip;

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



    bool PatchSignalGenerator::DoGenerate( Signal* aSignal )
    {

        InitializePatchArray();


        //n samples for event spacing.
        int PreEventCounter = 0;
        const int NPreEventSamples = 150000;

        std::thread Kassiopeia(KassiopeiaInit, gxml_filename);     // spawn new thread
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
                        DriveAntenna(PreEventCounter, index, aSignal);
                        PreEventCounter = 0; // reset
                    }
                    tLock.unlock();
                }

        }  // for loop

        printf("finished signal loop\n");

        fRunInProgress = false;  // tell Kassiopeia to finish.
        fDoneWithSignalGeneration = true;  // tell LMCCyclotronRadExtractor
        //if (fEventInProgress)
        //  if (ReceivedKassReady())
        WakeBeforeEvent();
        Kassiopeia.join();

        return true;
    }

} /* namespace locust */

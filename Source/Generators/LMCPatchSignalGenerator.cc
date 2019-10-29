/*
 * LMCPatchSignalGenerator.cc
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
 */

#include "LMCPatchSignalGenerator.hh"
#include "LMCRunKassiopeia.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>

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
		fZShiftArray( 0. ),
        fPatchSpacing( 0. ),
        gxml_filename("blank.xml"),
		fTextFileWriting( 0 ),
        fphiLO(0.),
		fNPreEventSamples( 150000 ),
        EFieldBuffer( 1 ),
        EPhaseBuffer( 1 ),
        EAmplitudeBuffer( 1 ),
        EFrequencyBuffer( 1 ),
        LOPhaseBuffer( 1 ),
        IndexBuffer( 1 ),
        PatchFIRBuffer( 1 ),
        fFieldBufferSize( 50 ),
		fInterface( new KassLocustInterface() )

    {

        fRequiredSignalState = Signal::kTime;

        KLInterfaceBootstrapper::get_instance()->SetInterface( fInterface );

    }

    PatchSignalGenerator::~PatchSignalGenerator()
    {
    }

    bool PatchSignalGenerator::Configure( const scarab::param_node& aParam )
    {
    	if(!fReceiverFIRHandler.Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver FIRHandler class");
    	}

    	if(!fPowerCombiner.Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver PowerCombiner class");
    	}

        if( aParam.has( "buffer-size" ) )
        {
        	fFieldBufferSize = aParam["buffer-size"]().as_int();
        	fHilbertTransform.SetBufferSize(aParam["buffer-size"]().as_int());
        }

    	if(!fHilbertTransform.Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver HilbertTransform class");
    	}

    	if( aParam.has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam["lo-frequency"]().as_double();
        }

        if( aParam.has( "array-radius" ) )
        {
            fArrayRadius = aParam["array-radius"]().as_double();
        }
        if( aParam.has( "npatches-per-strip" ) )
        {
            fNPatchesPerStrip = aParam["npatches-per-strip"]().as_int();
        }
        if( aParam.has( "patch-spacing" ) )
        {
            fPatchSpacing = aParam["patch-spacing"]().as_double();
        }
        if( aParam.has( "zshift-array" ) )
        {
            fZShiftArray = aParam["zshift-array"]().as_double();
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

    void PatchSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    void PatchSignalGenerator::KassiopeiaInit(const std::string &aFile)
    {
        RunKassiopeia tRunKassiopeia;
        tRunKassiopeia.Run(aFile, fInterface);
        return;
    }


    bool PatchSignalGenerator::WakeBeforeEvent()
    {
        fInterface->fPreEventCondition.notify_one();
        return true;
    }

    bool PatchSignalGenerator::ReceivedKassReady()
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        printf("LMC about to wait ..\n ");

        if((fInterface->fRunInProgress)&&( ! fInterface->fKassEventReady))
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


    double PatchSignalGenerator::GetAOIFactor(LMCThreeVector IncidentKVector, double PatchPhi)
    {
        LMCThreeVector PatchNormalVector;
        PatchNormalVector.SetComponents(cos(PatchPhi), sin(PatchPhi), 0.0);
        double AOIFactor = fabs(IncidentKVector.Unit().Dot(PatchNormalVector));
        //printf("cos aoi is %f\n", AOIFactor);
        return AOIFactor;
    }


    // fields incident on patch.
    void PatchSignalGenerator::RecordIncidentFields(FILE *fp, LMCThreeVector IncidentMagneticField, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double PatchPhi, double DopplerFrequency)
    {
        double AOIFactor = GetAOIFactor(IncidentKVector, PatchPhi);  // k dot patchnormal
        LMCThreeVector PatchPolarizationVector;
        PatchPolarizationVector.SetComponents(-sin(PatchPhi), cos(PatchPhi), 0.0);
        LMCThreeVector PatchCrossPolarizationVector;
        PatchCrossPolarizationVector.SetComponents(0.0, 0.0, 1.0);  // axial direction.

        double EFieldPatch = IncidentElectricField.Dot(PatchPolarizationVector);
        double BFieldPatch = IncidentMagneticField.Dot(PatchCrossPolarizationVector);
        fprintf(fp, "%10.4g %10.4g %10.4g %10.4g\n", EFieldPatch, BFieldPatch, DopplerFrequency/2./LMCConst::Pi(), fInterface->fTOld);
    }



    double PatchSignalGenerator::GetFIRSample(int nfilterbins, double dtfilter, unsigned channel, unsigned patch)
    {

    	double fieldfrequency = EFrequencyBuffer[channel*fNPatchesPerStrip+patch].front();
    	double HilbertMag = 0.;
    	double HilbertPhase = 0.;
    	double convolution = 0.0;

    	if (fabs(EFieldBuffer[channel*fNPatchesPerStrip+patch].front()) > 0.)  // field arrived yet?
    	{

    		double* HilbertMagPhaseMean = new double[3];
    		HilbertMagPhaseMean = fHilbertTransform.GetMagPhaseMean(EFieldBuffer[channel*fNPatchesPerStrip+patch], EFrequencyBuffer[channel*fNPatchesPerStrip+patch]);
    		HilbertMag = HilbertMagPhaseMean[0];
    		HilbertPhase = HilbertMagPhaseMean[1];
    		delete[] HilbertMagPhaseMean;

    		for (int i=0; i < nfilterbins; i++)  // populate filter with field.
    		{
    			HilbertPhase += 2.*3.1415926*fieldfrequency*dtfilter;
    			PatchFIRBuffer[channel*fNPatchesPerStrip+patch].push_back(HilbertMag*cos(HilbertPhase));
    			PatchFIRBuffer[channel*fNPatchesPerStrip+patch].pop_front();
    		}

    		convolution=fReceiverFIRHandler.ConvolveWithFIRFilter(PatchFIRBuffer[channel*fNPatchesPerStrip+patch]);

    		PatchFIRBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();  // memory deallocation.
    		return convolution;
    	}
    	else return 0.;

    }



    // EField cross pol with aoi dot product, at patch.
    double PatchSignalGenerator::GetEFieldCoPol(PatchAntenna* currentPatch, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double PatchPhi, double DopplerFrequency)
    {
        double AOIFactor = GetAOIFactor(IncidentKVector, PatchPhi);  // k dot patchnormal
        LMCThreeVector PatchPolarizationVector = currentPatch->GetPolarizationDirection();
        double EFieldCoPol = IncidentElectricField.Dot(PatchPolarizationVector) * AOIFactor;

        return EFieldCoPol;
    }



    void PatchSignalGenerator::DriveAntenna(FILE *fp, int PreEventCounter, unsigned index, Signal* aSignal, int nfilterbins, double dtfilter)
    {

        const int signalSize = aSignal->TimeSize();
        unsigned sampleIndex = 0;

        //Receiver Properties
        fphiLO += 2. * LMCConst::Pi() * fLO_Frequency * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        double tReceiverTime = fInterface->fTOld;

        PatchAntenna *currentPatch;
        unsigned tTotalPatchIndex = 0;

        fFieldSolver.SetParticleHistory(fInterface->fParticleHistory);

        for(int channelIndex = 0; channelIndex < allChannels.size(); ++channelIndex)
        {
            double PatchPhi = (double)channelIndex*360./allChannels.size()*LMCConst::Pi()/180.; // radians.    
            for(int patchIndex = 0; patchIndex < allChannels[channelIndex].size(); ++patchIndex)
            {
                currentPatch = &allChannels[channelIndex][patchIndex]; 
                sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample

                fFieldSolver.SetFieldEvent(tReceiverTime, tTotalPatchIndex);
                fFieldSolver.SolveFieldSolutions();

                LMCThreeVector tRadiatedElectricField = fFieldSolver.GetElectricField();
                LMCThreeVector tRadiatedMagneticField = fFieldSolver.GetMagneticField();
                locust::Particle tCurrentParticle = fFieldSolver.GetRetardedParticle();

                LMCThreeVector tDirection = currentPatch->GetPosition() - tCurrentParticle.GetPosition(true);
                double tVelZ = tCurrentParticle.GetVelocity(true).Z();
                double tCosTheta =  tVelZ * tDirection.Z() /  tDirection.Magnitude() / fabs(tVelZ);
                double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency() / ( 1. - fabs(tVelZ) / LMCConst::C() * tCosTheta);

 		        double tEFieldCoPol = GetEFieldCoPol(currentPatch, tRadiatedElectricField, tRadiatedElectricField.Cross(tRadiatedMagneticField), PatchPhi, tDopplerFrequency);
                if (fTextFileWriting==1) RecordIncidentFields(fp, tRadiatedMagneticField, tRadiatedElectricField, tRadiatedElectricField.Cross(tRadiatedMagneticField) , PatchPhi, tDopplerFrequency);

 	            FillBuffers(aSignal, tDopplerFrequency, tEFieldCoPol, fphiLO, index, channelIndex, patchIndex);
 	            double VoltageFIRSample = GetFIRSample(nfilterbins, dtfilter, channelIndex, patchIndex);
 	            fPowerCombiner.AddOneVoltageToStripSum(aSignal, VoltageFIRSample, fphiLO, patchIndex, IndexBuffer[channelIndex*fNPatchesPerStrip+patchIndex].front());
                PopBuffers(channelIndex, patchIndex);

                ++tTotalPatchIndex;

            } // patch loop

        } // channels loop


        fInterface->fTOld += 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());  // advance time here instead of in step modifier.  This preserves the freefield sampling.

    }

    void PatchSignalGenerator::FillBuffers(Signal* aSignal, double DopplerFrequency, double EFieldValue, double LOPhase, unsigned index, unsigned channel, unsigned patch)
    {
    	EFieldBuffer[channel*fNPatchesPerStrip+patch].push_back(EFieldValue);
    	EFrequencyBuffer[channel*fNPatchesPerStrip+patch].push_back(DopplerFrequency/2./LMCConst::Pi());
    	LOPhaseBuffer[channel*fNPatchesPerStrip+patch].push_back(LOPhase);
    	IndexBuffer[channel*fNPatchesPerStrip+patch].push_back(channel*aSignal->TimeSize()*aSignal->DecimationFactor() + index);
    }





    void PatchSignalGenerator::PopBuffers(unsigned channel, unsigned patch)
    {

    	EFieldBuffer[channel*fNPatchesPerStrip+patch].pop_front();
    	EFrequencyBuffer[channel*fNPatchesPerStrip+patch].pop_front();
    	LOPhaseBuffer[channel*fNPatchesPerStrip+patch].pop_front();
    	IndexBuffer[channel*fNPatchesPerStrip+patch].pop_front();
    	EFieldBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        EFrequencyBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        LOPhaseBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        IndexBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();

    }




    void PatchSignalGenerator::InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize)
    {

    	FieldBuffer aFieldBuffer;
    	EFieldBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
    	EFrequencyBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
    	LOPhaseBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
    	IndexBuffer = aFieldBuffer.InitializeUnsignedBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
    	PatchFIRBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, filterbuffersize);

    }


    void PatchSignalGenerator::CleanupBuffers()
    {
    	FieldBuffer aFieldBuffer;
    	EFieldBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    	EFrequencyBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    	LOPhaseBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    	IndexBuffer = aFieldBuffer.CleanupBuffer(IndexBuffer);

    }



    bool PatchSignalGenerator::InitializePowerCombining()
    {
    	fPowerCombiner.SetSMatrixParameters(fNPatchesPerStrip);
    	fPowerCombiner.SetVoltageDampingFactors(fNPatchesPerStrip);
    	return true;

    }


    bool PatchSignalGenerator::InitializePatchArray()
    {

        if(!fReceiverFIRHandler.ReadFIRFile())
        {
            return false;
        }


        const unsigned nChannels = fNChannels;
        const int nReceivers = fNPatchesPerStrip;

        const double patchSpacingZ = fPatchSpacing;
        const double patchRadius = fArrayRadius;
        double zPosition;
        double theta;
        const double dThetaArray = 2. * LMCConst::Pi() / nChannels; //Divide the circle into nChannels
        const double dRotateVoltages = 0.;  // set to zero to not rotate patch polarities.

        PatchAntenna modelPatch;

        allChannels.resize(nChannels);

        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            theta = channelIndex * dThetaArray;

            for(int receiverIndex = 0; receiverIndex < nReceivers; ++receiverIndex)
            {
                zPosition =  fZShiftArray + (receiverIndex - (nReceivers - 1.) /2.) * patchSpacingZ;

                modelPatch.SetCenterPosition({patchRadius * cos(theta) , patchRadius * sin(theta) , zPosition }); 
                modelPatch.SetPolarizationDirection({sin(theta), -cos(theta), 0.});
           
                modelPatch.SetNormalDirection({-cos(theta), -sin(theta), 0.}); //Say normals point inwards
                allChannels[channelIndex].AddReceiver(modelPatch);
                fFieldSolver.AddFieldPoint(modelPatch.GetPosition());
            }
        }
        return true;
    }



    bool PatchSignalGenerator::DoGenerate( Signal* aSignal )
    {

        FILE *fp = fopen("incidentfields.txt", "w");

        InitializePatchArray();
        InitializePowerCombining();

        //n samples for event spacing.
        int PreEventCounter = 0;
        fInterface->fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        std::thread tKassiopeia (&PatchSignalGenerator::KassiopeiaInit, this, gxml_filename);     // spawn new thread
        fInterface->fRunInProgress = true;

        int nfilterbins = fReceiverFIRHandler.GetFilterSize();
        double dtfilter = fReceiverFIRHandler.GetFilterResolution();
        unsigned nfieldbufferbins = fFieldBufferSize;
        InitializeBuffers(nfilterbins, nfieldbufferbins);

        for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        {
            if ((!fInterface->fEventInProgress) && (fInterface->fRunInProgress) && (!fInterface->fPreEventInProgress))
            {
            	if (ReceivedKassReady()) fInterface->fPreEventInProgress = true;
            	fInterface->fPreEventInProgress = true;
            	printf("LMC says it ReceivedKassReady(), fRunInProgress is %d\n", fInterface->fRunInProgress);

            }

            if ((fInterface->fPreEventInProgress)&&(fInterface->fRunInProgress))
            {
                PreEventCounter += 1;

                if (PreEventCounter > fNPreEventSamples) // finished pre-samples.  Start event.
                {
                    fInterface->fPreEventInProgress = false;  // reset.
                    fInterface->fEventInProgress = true;
                    printf("LMC about to WakeBeforeEvent()\n");
                    WakeBeforeEvent();  // trigger Kass event.
                }
            }

            if (fInterface->fEventInProgress)  // fEventInProgress
            {
                    std::unique_lock< std::mutex >tLock( fInterface->fMutexDigitizer, std::defer_lock );
                    tLock.lock();
                    fInterface->fDigitizerCondition.wait( tLock );
                    if (fInterface->fEventInProgress)
                    {
                        DriveAntenna(fp, PreEventCounter, index, aSignal, nfilterbins, dtfilter);
                        PreEventCounter = 0; // reset
                    }
                    tLock.unlock();
            }
        }  // for loop

        fInterface->fRunInProgress = false;  // tell Kassiopeia to finish.
        fInterface->fDoneWithSignalGeneration = true;  // tell LMCCyclotronRadExtractor
        fclose(fp);
        printf("finished signal loop.\n");
        WakeBeforeEvent();  // trigger one last Kass event if we are locked up.
        tKassiopeia.join();  // finish thread


        return true;
    }

} /* namespace locust */


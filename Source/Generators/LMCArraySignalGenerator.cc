/*
 * LMCArraySignalGenerator.cc
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
 */

#include "LMCArraySignalGenerator.hh"
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
    LOGGER( lmclog, "ArraySignalGenerator" );

    MT_REGISTER_GENERATOR(ArraySignalGenerator, "array-signal");

    ArraySignalGenerator::ArraySignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fLO_Frequency( 0.),
        fArrayRadius( 0. ),
        fNElementsPerStrip( 0. ),
		fZShiftArray( 0. ),
        fElementSpacing( 0. ),
        gxml_filename("blank.xml"),
		fTextFileWriting( 0 ),
        fphiLO(0.),
        EFieldBuffer( 1 ),
        EPhaseBuffer( 1 ),
        EAmplitudeBuffer( 1 ),
        EFrequencyBuffer( 1 ),
        LOPhaseBuffer( 1 ),
        IndexBuffer( 1 ),
        ElementFIRBuffer( 1 ),
        fFieldBufferSize( 50 ),
		fSwapFrequency( 1000 )

    {
        fRequiredSignalState = Signal::kTime;
    }

    ArraySignalGenerator::~ArraySignalGenerator()
    {
    }

    bool ArraySignalGenerator::Configure( const scarab::param_node& aParam )
    {
    	if(!fTFReceiverHandler.Configure(aParam))
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
    		LERROR(lmclog,"Error configuring buffer sizes in receiver HilbertTransform class");
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
        if( aParam.has( "element-spacing" ) )
        {
            fElementSpacing = aParam["element-spacing"]().as_double();
        }
        if( aParam.has( "zshift-array" ) )
        {
            fZShiftArray = aParam["zshift-array"]().as_double();
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

    void ArraySignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
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



    bool ArraySignalGenerator::WakeBeforeEvent()
    {
        fPreEventCondition.notify_one();
        return true;
    }

    bool ArraySignalGenerator::ReceivedKassReady()
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
    void ArraySignalGenerator::RecordIncidentFields(FILE *fp,  double t_old, int elementIndex, double zelement, double tEFieldCoPol)
    {
    	if (t_old > 0.5e-9)
         	{
         	fprintf(fp, "%d %g %g\n", elementIndex, zelement, tEFieldCoPol);
         	}
    }



    double ArraySignalGenerator::GetFIRSample(int nfilterbins, double dtfilter, unsigned channel, unsigned element)
    {

    	double fieldfrequency = EFrequencyBuffer[channel*fNElementsPerStrip+element].front();
    	double HilbertMag = 0.;
    	double HilbertPhase = 0.;
    	double convolution = 0.0;

    	if (fabs(EFieldBuffer[channel*fNElementsPerStrip+element].front()) > 0.)  // field arrived yet?
    	{

    		double* HilbertMagPhaseMean = new double[3];
    		HilbertMagPhaseMean = fHilbertTransform.GetMagPhaseMean(EFieldBuffer[channel*fNElementsPerStrip+element], EFrequencyBuffer[channel*fNElementsPerStrip+element]);
    		HilbertMag = HilbertMagPhaseMean[0];
    		HilbertPhase = HilbertMagPhaseMean[1];
    		delete[] HilbertMagPhaseMean;

    		for (int i=0; i < nfilterbins; i++)  // populate filter with field.
    		{
    			HilbertPhase += 2.*3.1415926*fieldfrequency*dtfilter;
    			ElementFIRBuffer[channel*fNElementsPerStrip+element].push_back(HilbertMag*cos(HilbertPhase));
    			ElementFIRBuffer[channel*fNElementsPerStrip+element].pop_front();
    		}

    		convolution=fTFReceiverHandler.ConvolveWithFIRFilter(ElementFIRBuffer[channel*fNElementsPerStrip+element]);

    		return convolution;
    	}
    	else return 0.;

    }




    void ArraySignalGenerator::DriveAntenna(FILE *fp, int PreEventCounter, unsigned index, Signal* aSignal, int nfilterbins, double dtfilter)
    {

        const int signalSize = aSignal->TimeSize();
        unsigned sampleIndex = 0;
        const unsigned nChannels = fNChannels;
        const int nReceivers = fNElementsPerStrip;


        //Receiver Properties
        fphiLO += 2. * LMCConst::Pi() * fLO_Frequency * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        double tReceiverTime = t_old;

        unsigned tTotalElementIndex = 0;

        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            double ElementPhi = (double)channelIndex*360./nChannels*LMCConst::Pi()/180.; // radians.
            for(int elementIndex = 0; elementIndex < nReceivers; ++elementIndex)
            {
            	Receiver* currentElement = allRxChannels[channelIndex][elementIndex];

                sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample

                fFieldSolver.SetFieldEvent(tReceiverTime, tTotalElementIndex);
                fFieldSolver.SolveFieldSolutions();

                LMCThreeVector tRadiatedElectricField = fFieldSolver.GetElectricField();
                LMCThreeVector tRadiatedMagneticField = fFieldSolver.GetMagneticField();
                locust::Particle tCurrentParticle = fFieldSolver.GetRetardedParticle();

                LMCThreeVector tDirection = currentElement->GetPosition() - tCurrentParticle.GetPosition(true);
                double tVelZ = tCurrentParticle.GetVelocity(true).Z();
                double tCosTheta =  tVelZ * tDirection.Z() /  tDirection.Magnitude() / fabs(tVelZ);
                double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency() / ( 1. - fabs(tVelZ) / LMCConst::C() * tCosTheta);


 		        double tEFieldCoPol = tRadiatedElectricField.Dot(currentElement->GetPolarizationDirection())*currentElement->GetPatternFactor(tRadiatedElectricField.Cross(tRadiatedMagneticField), *currentElement);
                if (fTextFileWriting==1) RecordIncidentFields(fp, t_old, elementIndex, currentElement->GetPosition().GetZ(), tEFieldCoPol);

 	            FillBuffers(aSignal, tDopplerFrequency, tEFieldCoPol, fphiLO, index, channelIndex, elementIndex);
 	            double VoltageFIRSample = GetFIRSample(nfilterbins, dtfilter, channelIndex, elementIndex);
 	            fPowerCombiner.AddOneVoltageToStripSum(aSignal, VoltageFIRSample, fphiLO, elementIndex, IndexBuffer[channelIndex*fNElementsPerStrip+elementIndex].front());
                PopBuffers(channelIndex, elementIndex);

                ++tTotalElementIndex;

            } // element loop

        } // channels loop

        t_old += 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        if ( index%fSwapFrequency == 0 ) CleanupBuffers();  // release memory

    }

    void ArraySignalGenerator::FillBuffers(Signal* aSignal, double DopplerFrequency, double EFieldValue, double LOPhase, unsigned index, unsigned channel, unsigned element)
    {
    	EFieldBuffer[channel*fNElementsPerStrip+element].push_back(EFieldValue);
    	EFrequencyBuffer[channel*fNElementsPerStrip+element].push_back(DopplerFrequency/2./LMCConst::Pi());
    	LOPhaseBuffer[channel*fNElementsPerStrip+element].push_back(LOPhase);
    	IndexBuffer[channel*fNElementsPerStrip+element].push_back(channel*aSignal->TimeSize()*aSignal->DecimationFactor() + index);
    }





    void ArraySignalGenerator::PopBuffers(unsigned channel, unsigned element)
    {

    	EFieldBuffer[channel*fNElementsPerStrip+element].pop_front();
    	EFrequencyBuffer[channel*fNElementsPerStrip+element].pop_front();
    	LOPhaseBuffer[channel*fNElementsPerStrip+element].pop_front();
    	IndexBuffer[channel*fNElementsPerStrip+element].pop_front();

    }




    void ArraySignalGenerator::InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize)
    {
    	FieldBuffer aFieldBuffer;
    	EFieldBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
    	EFrequencyBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
    	LOPhaseBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
    	IndexBuffer = aFieldBuffer.InitializeUnsignedBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
    	ElementFIRBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, filterbuffersize);
    }


    void ArraySignalGenerator::CleanupBuffers()
    {
    	FieldBuffer aFieldBuffer;
    	EFieldBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    	EFrequencyBuffer = aFieldBuffer.CleanupBuffer(EFrequencyBuffer);
    	LOPhaseBuffer = aFieldBuffer.CleanupBuffer(LOPhaseBuffer);
    	ElementFIRBuffer = aFieldBuffer.CleanupBuffer(ElementFIRBuffer);
    	IndexBuffer = aFieldBuffer.CleanupBuffer(IndexBuffer);
    }



    bool ArraySignalGenerator::InitializePowerCombining()
    {
    	fPowerCombiner.SetSMatrixParameters(fNElementsPerStrip);
    	if (!fPowerCombiner.SetVoltageDampingFactors(fNElementsPerStrip, fElementSpacing) )
    	{
    		return false;
    	}
    	else
    	{
    		return true;
    	}
    }


    bool ArraySignalGenerator::InitializeElementArray()
    {

        if(!fTFReceiverHandler.ReadHFSSFile())
        {
            return false;
        }

        const unsigned nChannels = fNChannels;
        const int nReceivers = fNElementsPerStrip;

        const double elementSpacingZ = fElementSpacing;
        const double elementRadius = fArrayRadius;
        double zPosition;
        double theta;
        const double dThetaArray = 2. * LMCConst::Pi() / nChannels; //Divide the circle into nChannels
        const double dRotateVoltages = 0.;  // set to zero to not rotate element polarities.

        Receiver* modelElement = new Receiver;
        modelElement = fPowerCombiner.ChooseElement();

        allRxChannels.resize(nChannels);

        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            theta = channelIndex * dThetaArray;

            for(int receiverIndex = 0; receiverIndex < nReceivers; ++receiverIndex)
            {
                zPosition =  fZShiftArray + (receiverIndex - (nReceivers - 1.) /2.) * elementSpacingZ;

                if (fPowerCombiner.GetPowerCombiner() == 7)  // single patch
                {
                	zPosition = 0.;
                }

                	modelElement->SetCenterPosition({elementRadius * cos(theta) , elementRadius * sin(theta) , zPosition });
                	modelElement->SetPolarizationDirection({sin(theta), -cos(theta), 0.});
                	modelElement->SetNormalDirection({-cos(theta), -sin(theta), 0.}); //Say normals point inwards
                	allRxChannels[channelIndex].AddReceiver(modelElement);
                	fFieldSolver.AddFieldPoint(modelElement->GetPosition());
            }
        }

        return true;
    }



    bool ArraySignalGenerator::DoGenerate( Signal* aSignal )
    {

        FILE *fp = fopen("incidentfields.txt", "w");

        if(!InitializeElementArray())
        {
        	LERROR(lmclog,"Error configuring Element array");
            exit(-1);
        }


        if (!InitializePowerCombining() )
        {
        	LERROR(lmclog,"Error configuring Power Combining");
            exit(-1);
        }

        //n samples for event spacing.
        int PreEventCounter = 0;
        const int NPreEventSamples = 150000;
        fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        std::thread Kassiopeia(KassiopeiaInit, gxml_filename);     // spawn new thread
        fRunInProgress = true;

        int nfilterbins = fTFReceiverHandler.GetFilterSize();
        double dtfilter = fTFReceiverHandler.GetFilterResolution();
        unsigned nfieldbufferbins = fFieldBufferSize;
        InitializeBuffers(nfilterbins, nfieldbufferbins);

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
                        DriveAntenna(fp, PreEventCounter, index, aSignal, nfilterbins, dtfilter);
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

} /* namespace locust */


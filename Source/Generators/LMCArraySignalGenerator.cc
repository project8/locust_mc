/*
 * LMCArraySignalGenerator.cc
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
 */

#include "LMCArraySignalGenerator.hh"

#include "LMCRunKassiopeia.hh"

#include "logger.hh"

#include <chrono>
#include <thread>


namespace locust
{
    LOGGER( lmclog, "ArraySignalGenerator" );

    MT_REGISTER_GENERATOR(ArraySignalGenerator, "array-signal");

    ArraySignalGenerator::ArraySignalGenerator( const std::string& aName ) :
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
        EFieldBuffer( 1 ),
        EPhaseBuffer( 1 ),
        EAmplitudeBuffer( 1 ),
        EFrequencyBuffer( 1 ),
        LOPhaseBuffer( 1 ),
        IndexBuffer( 1 ),
        ElementFIRBuffer( 1 ),
        FIRfrequencyBuffer( 1 ),
        fFieldBufferSize( 50 ),
		fFIRzeroBuffer( 20 ),
		fSwapFrequency( 1000 ),
		fKassNeverStarted( false ),
		fSkippedSamples( false ),
		fAllowFastSampling( false ),
		fInterface( new KassLocustInterface() )
    {
        fRequiredSignalState = Signal::kTime;

        KLInterfaceBootstrapper::get_instance()->SetInterface( fInterface );
    }

    ArraySignalGenerator::~ArraySignalGenerator()
    {
    }

    bool ArraySignalGenerator::Configure( const scarab::param_node& aParam )
    {

    	if (aParam.has( "power-combining-feed" ))
    	{
    		int npowercombiners = 0;

        	if(aParam["power-combining-feed"]().as_string() == "voltage-divider")
        	{
        		npowercombiners += 1;
        		fPowerCombiner = new VoltageDivider;
        		if(!fPowerCombiner->Configure(aParam))
        		{
        			LERROR(lmclog,"Error configuring voltage divider.");
        			exit(-1);
        		}
        		fAntennaElementPositioner = new AntennaElementPositioner;
           		if(!fAntennaElementPositioner->Configure(aParam))
            	{
            		LERROR(lmclog,"Error configuring antenna element positioner.");
            		exit(-1);
            	}
        	}

        	if(aParam["power-combining-feed"]().as_string() == "slotted-waveguide")
        	{
        		npowercombiners += 1;
        		fPowerCombiner = new SlottedWaveguide;
        		if(!fPowerCombiner->Configure(aParam))
        		{
        			LERROR(lmclog,"Error configuring slotted waveguide.");
        			exit(-1);
        		}
        		fAntennaElementPositioner = new AntennaElementPositioner;
           		if(!fAntennaElementPositioner->Configure(aParam))
            	{
            		LERROR(lmclog,"Error configuring antenna element positioner.");
            		exit(-1);
            	}
        	}

        	if(aParam["power-combining-feed"]().as_string() == "single-patch")
        	{
        		npowercombiners += 1;
        		fPowerCombiner = new SinglePatch;
        		if(!fPowerCombiner->Configure(aParam))
        		{
        			LERROR(lmclog,"Error configuring single patch.");
        			exit(-1);
        		}
        		fAntennaElementPositioner = new SinglePatchPositioner;
           		if(!fAntennaElementPositioner->Configure(aParam))
            	{
            		LERROR(lmclog,"Error configuring single patch positioner.");
            		exit(-1);
            	}
        	}

        	if(aParam["power-combining-feed"]().as_string() == "corporate")
        	{
        		npowercombiners += 1;
        		fPowerCombiner = new CorporateFeed;
        		if(!fPowerCombiner->Configure(aParam))
        		{
        			LERROR(lmclog,"Error configuring corporate feed.");
        			exit(-1);
        		}
        		fAntennaElementPositioner = new AntennaElementPositioner;
           		if(!fAntennaElementPositioner->Configure(aParam))
            	{
            		LERROR(lmclog,"Error configuring antenna element positioner.");
            		exit(-1);
            	}
        	}

        	if(aParam["power-combining-feed"]().as_string() == "s-matrix")
        	{
        		npowercombiners += 1;
        		fPowerCombiner = new SMatrix;
        		if(!fPowerCombiner->Configure(aParam))
        		{
        			LERROR(lmclog,"Error configuring s matrix.");
        			exit(-1);
        		}
        		fAntennaElementPositioner = new AntennaElementPositioner;
           		if(!fAntennaElementPositioner->Configure(aParam))
            	{
            		LERROR(lmclog,"Error configuring antenna element positioner.");
            		exit(-1);
            	}
        	}


        	if((aParam["power-combining-feed"]().as_string() == "unit-cell-one-quarter")||
               (aParam["power-combining-feed"]().as_string() == "unit-cell-seven-eighths")||
               (aParam["power-combining-feed"]().as_string() == "unit-cell-nine-sixteenths"))
        	{
        		npowercombiners += 1;
        		fPowerCombiner = new UnitCell;
        		if(!fPowerCombiner->Configure(aParam))
        		{
        			LERROR(lmclog,"Error configuring unit cell.");
        			exit(-1);
        		}
        		fAntennaElementPositioner = new AntennaElementPositioner;
           		if(!fAntennaElementPositioner->Configure(aParam))
            	{
            		LERROR(lmclog,"Error configuring antenna element positioner.");
            		exit(-1);
            	}
        	}


        	if(aParam["power-combining-feed"]().as_string() == "series-feed")
        	{
        		npowercombiners += 1;
        		fPowerCombiner = new SeriesFeed;
        		if(!fPowerCombiner->Configure(aParam))
        		{
        			LERROR(lmclog,"Error configuring series feed.");
        		}
        		fAntennaElementPositioner = new AntennaElementPositioner;
           		if(!fAntennaElementPositioner->Configure(aParam))
            	{
            		LERROR(lmclog,"Error configuring antenna element positioner.");
            		exit(-1);
            	}
        	}


        	if (npowercombiners != 1)
        	{
        		LERROR(lmclog,"LMCArraySignalGenerator needs a single power combiner.  Please choose one value for power-combining-feed in the config file.");
                exit(-1);
        	}

    	}
        else
        {
    		LERROR(lmclog,"LMCArraySignalGenerator has been configured without a power combiner.  Please choose a value for power-combiner-feed in the config file.");
            exit(-1);
        }


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
        		LERROR(lmclog,"LMCArraySignalGenerator needs a single transmitter.  Please choose transmitter:antenna or transmitter:planewave or transmitter:kassiopeia in the config file.");
                exit(-1);
        	}
        }
        else
        {
    		LERROR(lmclog,"LMCArraySignalGenerator has been configured without a transmitter.  Please choose transmitter:antenna or transmitter:planewave or transmitter:kassiopeia in the config file.");
            exit(-1);
        }


    	if(!fTFReceiverHandler.Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver FIRHandler class");
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
        if( aParam.has( "allow-fast-sampling" ) )
        {
            fAllowFastSampling = aParam["allow-fast-sampling"]().as_bool();
        }

        return true;
    }

    void ArraySignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }



    void ArraySignalGenerator::KassiopeiaInit(const std::string &aFile)
    {
        RunKassiopeia tRunKassiopeia;
        tRunKassiopeia.Run(aFile, fInterface);
        return;
    }

    void ArraySignalGenerator::InitializeFieldPoints(std::vector< Channel<Receiver*> > allRxChannels)
    {
	for(int channelIndex = 0; channelIndex < fNChannels; ++channelIndex)
	{
            for(int elementIndex = 0; elementIndex < fNElementsPerStrip; ++elementIndex)
            {
            	fTransmitter->InitializeFieldPoint(allRxChannels[channelIndex][elementIndex]->GetPosition());
            }
	}
    }

    void ArraySignalGenerator::WakeBeforeEvent()
    {
        fInterface->fPreEventCondition.notify_one();
        return;
    }

    bool ArraySignalGenerator::ReceivedKassReady()
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


    // fields incident on element.
    void ArraySignalGenerator::RecordIncidentFields(FILE *fp,  double t_old, int elementIndex, double zelement, double tEFieldCoPol)
    {
    	if (t_old > 0.5e-9)
         	{
         	fprintf(fp, "%d %g %g\n", elementIndex, zelement, tEFieldCoPol);
         	}
    }



    double ArraySignalGenerator::GetFIRSample(int nFilterBinsRequired, double dtFilter, unsigned channel, unsigned element)
    {

    	double fieldfrequency = EFrequencyBuffer[channel*fNElementsPerStrip+element].front();
    	double HilbertMag = 0.;
    	double HilbertPhase = 0.;
    	double convolution = 0.0;

    	if (EFieldBuffer[channel*fNElementsPerStrip+element].front() != 0.)  // field arrived yet?
    	{

    		std::vector<double> HilbertMagPhaseMean; HilbertMagPhaseMean.resize(3);
    		HilbertMagPhaseMean = fHilbertTransform.GetMagPhaseMean(EFieldBuffer[channel*fNElementsPerStrip+element], EFrequencyBuffer[channel*fNElementsPerStrip+element]);
    		HilbertMag = HilbertMagPhaseMean[0];
    		HilbertPhase = HilbertMagPhaseMean[1];

            ////////////////////////////
            //text file for testing.   
            std::ofstream hilbertfile;
            hilbertfile.open("hilbertfile.txt", std::ios::app);
            hilbertfile << HilbertMag << ", " << HilbertPhase;
            hilbertfile << "\n";
            hilbertfile.close();
            ////////////////////////////

            ////////////////////////////
            //text file for testing.   
            std::ofstream fieldfile;
            fieldfile.open("fieldfile.txt", std::ios::app);
            fieldfile << EFieldBuffer[channel*fNElementsPerStrip+element].at(fFieldBufferSize/2);
            fieldfile << "\n";
            fieldfile.close();
            ////////////////////////////

    		// populate FIR filter with frequency for just this sample interval:
    		for (int i=0; i < nFilterBinsRequired; i++)
    		{
    			FIRfrequencyBuffer[channel*fNElementsPerStrip+element].push_back(fieldfrequency);
    			FIRfrequencyBuffer[channel*fNElementsPerStrip+element].pop_front();
    		}

    		// populate entire FIR filter with field, using frequencies from recent previous samples:
    		std::deque<double>::iterator it = FIRfrequencyBuffer[channel*fNElementsPerStrip+element].begin();
    		while (it != FIRfrequencyBuffer[channel*fNElementsPerStrip+element].end())
    		{
    			HilbertPhase += 2.*3.1415926*(*it)*dtFilter;
    			if (*it != 0.)
    			{
    				ElementFIRBuffer[channel*fNElementsPerStrip+element].push_back(HilbertMag*cos(HilbertPhase));
    			}
    			else
    			{
    				ElementFIRBuffer[channel*fNElementsPerStrip+element].push_back(0.);
    			}
    			ElementFIRBuffer[channel*fNElementsPerStrip+element].pop_front();
    			*it++;
    		}

    		convolution=fTFReceiverHandler.ConvolveWithFIRFilter(ElementFIRBuffer[channel*fNElementsPerStrip+element]);
    		return convolution;

    	}
    	else return 0.;

    }


    bool ArraySignalGenerator::DriveAntenna(FILE *fp, int startingIndex, unsigned index, Signal* aSignal, int nFilterBinsRequired, double dtFilter)
    {

        const int signalSize = aSignal->TimeSize();
        unsigned sampleIndex = 0;
        const unsigned nChannels = fNChannels;
        const int nReceivers = fNElementsPerStrip;

        //Receiver Properties
        fphiLO += 2. * LMCConst::Pi() * fLO_Frequency * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        double tReceiverTime =  fInterface->fTOld;

        unsigned tTotalElementIndex = 0;

        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            double ElementPhi = (double)channelIndex*360./nChannels*LMCConst::Pi()/180.; // radians.
            for(int elementIndex = 0; elementIndex < nReceivers; ++elementIndex)
            {
            	Receiver* currentElement = allRxChannels[channelIndex][elementIndex];
                sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample

                std::vector<double> tFieldSolution; //tFieldSolution.resize(2);
                if (!fTransmitter->IsKassiopeia())
                {
                	tFieldSolution = fTransmitter->GetEFieldCoPol(tTotalElementIndex, 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor()));
                }
                else
                {
                	tFieldSolution = fTransmitter->SolveKassFields(currentElement->GetPosition(), currentElement->GetPolarizationDirection(), tReceiverTime, tTotalElementIndex);
                }

                tFieldSolution[0] *= currentElement->GetPatternFactor(fTransmitter->GetIncidentKVector(tTotalElementIndex), *currentElement);

                if (fTextFileWriting==1) RecordIncidentFields(fp,  fInterface->fTOld, elementIndex, currentElement->GetPosition().GetZ(), tFieldSolution[1]);

 	            FillBuffers(aSignal, tFieldSolution[1], tFieldSolution[0], fphiLO, index, channelIndex, elementIndex);
 	            double VoltageFIRSample = GetFIRSample(nFilterBinsRequired, dtFilter, channelIndex, elementIndex);
            	if ((VoltageFIRSample == 0.)&&(index-startingIndex > fFieldBufferSize + fFIRzeroBuffer))
            	{
                    LERROR(lmclog,"A digitizer sample was skipped due to likely unresponsive thread.\n");
            		return false;
            	}

            	fPowerCombiner->AddOneVoltageToStripSum(aSignal, VoltageFIRSample, fphiLO, elementIndex, IndexBuffer[channelIndex*fNElementsPerStrip+elementIndex].front());
                PopBuffers(channelIndex, elementIndex);

                ++tTotalElementIndex;

            } // element loop

        } // channels loop

        fInterface->fTOld += 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        if ( index%fSwapFrequency == 0 ) CleanupBuffers();  // release memory
        return true;

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
    	FIRfrequencyBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, filterbuffersize);
    }


    void ArraySignalGenerator::CleanupBuffers()
    {
    	FieldBuffer aFieldBuffer;
    	EFieldBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    	EFrequencyBuffer = aFieldBuffer.CleanupBuffer(EFrequencyBuffer);
    	LOPhaseBuffer = aFieldBuffer.CleanupBuffer(LOPhaseBuffer);
    	ElementFIRBuffer = aFieldBuffer.CleanupBuffer(ElementFIRBuffer);
    	FIRfrequencyBuffer = aFieldBuffer.CleanupBuffer(FIRfrequencyBuffer);
    	IndexBuffer = aFieldBuffer.CleanupBuffer(IndexBuffer);
    }


    bool ArraySignalGenerator::InitializeElementArray()
    {

        if(!fTFReceiverHandler.ReadHFSSFile())
        {
            return false;
        }

        if(!fHilbertTransform.SetupHilbertTransform())
        {
            LERROR(lmclog,"Error initializing Hilbert transform.\n");
            return false;
        }

        const unsigned nChannels = fNChannels;
        const unsigned nSubarrays = fNSubarrays;
        const int nReceivers = fNElementsPerStrip;

        const double elementSpacingZ = fElementSpacing;
        const double elementRadius = fArrayRadius;
        double zPosition;
        double theta;
        const double dThetaArray = 2. * LMCConst::Pi() / (nChannels/nSubarrays); //Divide the circle into nChannels
        const double dRotateVoltages = 0.;  // set to zero to not rotate element polarities.

        allRxChannels.resize(nChannels);

        	for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        	{
        		theta = fAntennaElementPositioner->GetTheta(channelIndex, dThetaArray);

        		for(int receiverIndex = 0; receiverIndex < nReceivers; ++receiverIndex)
        		{
        			zPosition = fAntennaElementPositioner->GetPositionZ(fZShiftArray, channelIndex, nChannels,
        						nSubarrays, nReceivers, elementSpacingZ, receiverIndex);

        			Receiver* modelElement = fPowerCombiner->ChooseElement();  // patch or slot?

        			fAntennaElementPositioner->PlaceElement(*modelElement, elementRadius, theta, zPosition);

        			allRxChannels[channelIndex].AddReceiver(modelElement);
        		}
        	}

        return true;
    }



    bool ArraySignalGenerator::DoGenerate( Signal* aSignal )
    {

        if(!InitializeElementArray())
        {
        	LERROR(lmclog,"Error configuring Element array");
            exit(-1);
        }

        FILE *fp;
        if (fTextFileWriting==1) fp = fopen("incidentfields.txt", "w");


        //n samples for event spacing in Kass.
        int PreEventCounter = 0;

        int nFilterBins = fTFReceiverHandler.GetFilterSize();
        double dtFilter = fTFReceiverHandler.GetFilterResolution();
        int nFilterBinsRequired = std::min( 1. / (fAcquisitionRate*1.e6*aSignal->DecimationFactor()) / dtFilter, (double)nFilterBins );
        if (!fAllowFastSampling) nFilterBinsRequired = nFilterBins;
        unsigned nFieldBufferBins = fFieldBufferSize;
        InitializeBuffers(nFilterBins, nFieldBufferBins);

        InitializeFieldPoints(allRxChannels);

        if (!fTransmitter->IsKassiopeia())
        {
        	for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        	{
        		DriveAntenna(fp, PreEventCounter, index, aSignal, nFilterBinsRequired, dtFilter);
        	}  // for loop
        	return true;
        }


        if (fTransmitter->IsKassiopeia())
        {
        	bool fTruth = false;
        	int startingIndex;
            fInterface->fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        	std::thread tKassiopeia (&ArraySignalGenerator::KassiopeiaInit, this, gxml_filename); // spawn new thread

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
                    		if (DriveAntenna(fp, startingIndex, index, aSignal, nFilterBinsRequired, dtFilter))
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


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
        fDoGenerateFunc( &CavitySignalGenerator::DoGenerateTime ),
        fLO_Frequency( 0. ),
		fDeltaT( 0. ),
        gxml_filename("blank.xml"),
        fphiLO(0.),
        fNPreEventSamples( 150000 ),
        fRandomPreEventSamples( false ),
		fTrackDelaySeed( 0 ),
        fThreadCheckTime(100),
        fKassNeverStarted( false ),
        fAliasedFrequencies( false ),
        fOverrideAliasing( false ),
        fBypassTF( false ),
        fNormCheck( false ),
        fIntermediateFile( false ),
        fUseDirectKassPower( true ),
        fAliasingIsChecked( false ),
		fUnitTestRootFile( false ),
        fInterface( nullptr )  // Initialize fInterface to (nullptr) instead of to (new KassLocustInterface())
    {
    }

    CavitySignalGenerator::~CavitySignalGenerator()
    {
    }


    bool CavitySignalGenerator::ConfigureInterface( Signal* aSignal )
    {
    	if ( fInterface == nullptr ) fInterface.reset( new KassLocustInterface() );
        KLInterfaceBootstrapper::get_instance()->SetInterface( fInterface );

    	const scarab::param_node& tParam = *GetParameters();

    	if( tParam.has( "rectangular-waveguide" ) )
        {
        	fInterface->fbWaveguide = tParam["rectangular-waveguide"]().as_bool();
        	if ( tParam.has ( "direct-kass-power" ) )
        	{
        		fUseDirectKassPower = tParam["direct-kass-power"]().as_bool();
        	}
        }

        // Select mode fields and power combiners:
        if (fInterface->fbWaveguide) // Waveguide
        {
        	fInterface->fField = new RectangularWaveguide;
        	fPowerCombiner = new WaveguideModes;
        }
        else // Cavity
        {
            if( tParam.has( "rectangular-cavity" ) )
            {
            	if ( tParam["rectangular-cavity"]().as_bool() ) fInterface->fField = new RectangularCavity;
            }
            else fInterface->fField = new CylindricalCavity;
        	fPowerCombiner = new CavityModes;
        }

        // Configure selected mode fields and power combiners:
        if (!fInterface->fField->Configure(tParam))
        {
        	LERROR(lmclog,"Error configuring LMCField.");
        	exit(-1);
        	return false;
        }
        if(!fPowerCombiner->Configure(tParam))
		{
			LERROR(lmclog,"Error configuring PowerCombiner.");
			exit(-1);
			return false;
		}
        fModeSet = fInterface->fField->ModeSelect(fInterface->fbWaveguide, 0);


        // Define generic response function for use in cavity or waveguide:
    	fTFReceiverHandler = new TFReceiverHandler;
    	if(!fTFReceiverHandler->Configure(tParam))
    	{
    		LERROR(lmclog,"Error configuring receiver FIRHandler class");
    		exit(-1);
    		return false;
    	}

    	// Configure the generic response function:
        if( tParam.has( "tf-receiver-filename" ) && tParam.has( "tf-receiver-bin-width" ) ) // If using HFSS output for cavity or waveguide
        {
        	if (( fUseDirectKassPower ) && ( fInterface->fbWaveguide == true ))
        	{
        		LWARN(lmclog,"Using direct Kass energy budget, and not HFSS data, due to parameter \"direct-kass-power\" = true");
        	}
        	else if (!fTFReceiverHandler->ReadHFSSFile(1,0,1,0)) // Placeholder values for mode indices.
        	{
        		LERROR(lmclog,"FIR has not been generated.");
        		exit(-1);
        		return false;
        	}
        }
        else if (!fInterface->fbWaveguide) // Cavity config follows
        {
        	// No HFSS file is present:  Cavity as damped harmonic oscillator

        	fAnalyticResponseFunction = new DampedHarmonicOscillator();
    		if ( !fAnalyticResponseFunction->Configure(tParam) || !(CrossCheckCavityConfig()) )
    		{
    			LERROR(lmclog,"DampedHarmonicOscillator was not configured.");
    			exit(-1);
    			return false;
    		}
    		if (!fTFReceiverHandler->ConvertAnalyticGFtoFIR(fModeSet, fAnalyticResponseFunction->GetGFarray(fModeSet)))
    		{
    			LERROR(lmclog,"GF->FIR was not generated.");
    			exit(-1);
    			return false;
    		}
        } // tParam.has( "tf-receiver-filename" )
        else
        {
        	// Can't find the required config information:
        	LERROR(lmclog, "For the rectangular waveguide, either \"direct-kass-power\" should be true, or both a \"tf-receiver-filename\" and \"tf-receiver-bin-width\" should be specified.");
        	exit(-1);
        	return false;
        }



    	// Configure back reaction:
    	if (( fInterface->fbWaveguide ) && ( fPowerCombiner->GetWaveguideShortIsPresent() ))
    	{
    		if ( tParam.has( "back-reaction" ) )
    		{
    			// optional to switch off:
    			fInterface->fBackReaction = tParam["back-reaction"]().as_bool();
        		LWARN(lmclog,"Setting back-reaction to " << fInterface->fBackReaction );
    		}
    		else
    		{
    			// if the short is present in the waveguide, default is for back reaction = true
    			fInterface->fBackReaction = true;
    		}
    	}
    	else if ( !fInterface->fbWaveguide ) // Cavity
    	{
    		if ( tParam.has( "back-reaction" ) )
    		{
    			// optional to switch off
    			fInterface->fBackReaction = tParam["back-reaction"]().as_bool();
        		LWARN(lmclog,"Setting back-reaction to " << fInterface->fBackReaction );
    		}
    		else
    		{
    			// default is true in the cavity
    			fInterface->fBackReaction = true;
    		}
    	}
		else
		{
			// No short, no cavity -> no back reaction.
			fInterface->fBackReaction = false;
    		LWARN(lmclog,"Setting back-reaction to " << fInterface->fBackReaction );
		}

        // Configure the transmitter, which defines the impulses taken from Kassiopeia.
        if( tParam.has( "transmitter" ))
        {
        	int ntransmitters = 0;

        	if(tParam["transmitter"]().as_string() == "kass-current")
        	{
        		ntransmitters += 1;
        		fInterface->fTransmitter = new KassCurrentTransmitter;
        		if(!fInterface->fTransmitter->Configure(tParam))
        		{
        			LERROR(lmclog,"Error Configuring kass-current transmitter class");
        		}
        	}

        	if (ntransmitters != 1)
        	{
        		LERROR(lmclog,"LMCCavitySignalGenerator needs the kass-current transmitter.  Please choose transmitter:kass-current in the config file.");
                exit(-1);
        	}
        }
        else
        {
    		LERROR(lmclog,"LMCCavitySignalGenerator needs the kass-current transmitter.  Please choose transmitter:kass-current in the config file.");
            exit(-1);
        }

        if( tParam.has( "trigger-confirm" ) )
        {
        	fInterface->fTriggerConfirm = tParam["trigger-confirm"]().as_int();
        }
        fInterface->fFastRecordLength = fRecordLength * aSignal->DecimationFactor();

        // Configure Locust-Kass interface classes and parameters:
        fInterface->fConfigureKass = new ConfigureKass();
        fInterface->fConfigureKass->SetParameters( tParam );

        fFieldCalculator = new FieldCalculator();
        if(!fFieldCalculator->Configure(tParam))
        {
            LERROR(lmclog,"Error configuring receiver FieldCalculator class from CavitySignal.");
        }
        fFieldCalculator->SetNFilterBinsRequired( 1. / (fAcquisitionRate*1.e6*aSignal->DecimationFactor()) );
        for (int mu=0; mu < fModeSet.size(); mu++)
        {
            fFieldCalculator->SetFilterSize( fTFReceiverHandler->GetFilterSizeArray(fModeSet[mu][0], fModeSet[mu][1], fModeSet[mu][2], fModeSet[mu][3]));
        }

    	return true;
    }


    void CavitySignalGenerator::SetParameters( const scarab::param_node& aParam )
    {
    	fParam = &aParam;
    }


    const scarab::param_node* CavitySignalGenerator::GetParameters()
    {
    	return fParam;
    }


    bool CavitySignalGenerator::Configure( const scarab::param_node& aParam )
    {

    	SetParameters( aParam );

    	// Configure signal parameters:
        if( aParam.has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam["lo-frequency"]().as_double();
        }
        if( aParam.has( "event-spacing-samples" ) )
        {
            fNPreEventSamples = aParam["event-spacing-samples"]().as_int();
            if (aParam.has( "random-spacing-samples" ))
            {
                if (aParam["random-spacing-samples"]().as_bool() == true)
                {
                    fRandomPreEventSamples = true;
                }
        	    if ( aParam.has( "random-track-seed" ) )
        	    {
        	    	// Randomize seed for start time using seed for track length as input:
        	        scarab::param_node default_setting;
        	        default_setting.add("name","uniform");
        	        std::shared_ptr< BaseDistribution> tSeedDistribution = fDistributionInterface.get_dist(default_setting);
        	        fDistributionInterface.SetSeed( aParam["random-track-seed"]().as_int() );
        	        int tSeed = 1.e8 * tSeedDistribution->Generate();
        	        SetSeed( tSeed );
        	    }
        	    else
        	    {
        	    	// Delay SetSeed to allow time stamp to advance between randomized tracks.
            		std::this_thread::sleep_for(std::chrono::milliseconds(2000));
        	    	SetSeed (time(NULL) );
        	    }
            }
        }
        if( aParam.has( "override-aliasing" ) )
        {
            fOverrideAliasing = aParam["override-aliasing"]().as_bool();
        }
        if( aParam.has( "bypass-tf" ) )
        {
        	fBypassTF = aParam["bypass-tf"]().as_bool();
        }
        if( aParam.has( "norm-check" ) )
        {
        	fNormCheck = aParam["norm-check"]().as_bool();
        }
        if( aParam.has( "intermediate-file" ) )
        {
        	fIntermediateFile = aParam["intermediate-file"]().as_bool();
        }
        if( aParam.has( "unit-test-root-file" ) )
        {
        	fUnitTestRootFile = aParam["unit-test-root-file"]().as_bool();
        }
        if( aParam.has( "xml-filename" ) )
        {
            gxml_filename = aParam["xml-filename"]().as_string();
        }

        return true;
    }

    bool CavitySignalGenerator::RecordRunParameters( Signal* aSignal )
    {
#ifdef ROOT_FOUND
    	fInterface->aRunParameter = new RunParameters();
    	fInterface->aRunParameter->fSamplingRateMHz = fAcquisitionRate;
    	fInterface->aRunParameter->fDecimationFactor = aSignal->DecimationFactor();
    	fInterface->aRunParameter->fLOfrequency = fLO_Frequency;
    	fInterface->aRunParameter->fRandomSeed = fTrackDelaySeed;
#endif
    	return true;
    }



    bool CavitySignalGenerator::SetSeed(int aSeed)
    {
        LPROG(lmclog,"Setting random seed for track delay to " << aSeed);
        fTrackDelaySeed = aSeed;
        return true;
    }


    bool CavitySignalGenerator::RandomizeStartDelay()
    {
        scarab::param_node default_setting;
        default_setting.add("name","uniform");
        fStartDelayDistribution = fDistributionInterface.get_dist(default_setting);
        fDistributionInterface.SetSeed( fTrackDelaySeed );
        int tNPreEventSamples = fNPreEventSamples * fStartDelayDistribution->Generate();
        LPROG(lmclog,"Randomizing the start delay to " << tNPreEventSamples << " fast samples.");
        fNPreEventSamples = tNPreEventSamples;

        return true;
    }

    bool CavitySignalGenerator::CrossCheckAliasing( Signal* aSignal, double dopplerFrequency )
    {
    	LPROG(lmclog,"Running some cross-checks ...");

    	AliasingUtility anAliasingUtility;
    	double fs = fAcquisitionRate;
    	double dr = aSignal->DecimationFactor();

    	if (!anAliasingUtility.CheckAliasing( dopplerFrequency, fLO_Frequency, fs, dr ))
    	{
    		return false;
    	}
    	else
    	{
    		return true;
    	}
    }

    bool CavitySignalGenerator::CrossCheckCavityConfig()
	{

    	LPROG(lmclog,"Running some cavity cross-checks ...");

        for (int mu=0; mu<fModeSet.size(); mu++)
        {
            bool bTE = fModeSet[mu][0];
            int l = fModeSet[mu][1];
            int m = fModeSet[mu][2];
            int n = fModeSet[mu][3];
            CavityUtility aCavityUtility;
            double timeResolution = fAnalyticResponseFunction->GetDHOTimeResolution(bTE, l, m, n);
            double thresholdFactor = fAnalyticResponseFunction->GetDHOThresholdFactor(bTE, l, m, n);
            double cavityFrequency = fAnalyticResponseFunction->GetCavityFrequency(bTE, l, m, n);
            double qExpected = fAnalyticResponseFunction->GetCavityQ(bTE, l, m, n);
            aCavityUtility.SetOutputFile(fUnitTestRootFile);
            if (!aCavityUtility.CheckCavityQ( bTE, l, m, n, timeResolution, thresholdFactor, cavityFrequency, qExpected ))
            {
                LERROR(lmclog,"The cavity Q does not look quite right.  Please tune the configuration "
                		"with the unit test as in bin/testLMCCavity [-h]");
                return false;
            }
        }

        return true;
	}

    void CavitySignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    Signal::State CavitySignalGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void CavitySignalGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &CavitySignalGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &CavitySignalGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }




    bool CavitySignalGenerator::DriveMode(Signal* aSignal, unsigned index)
    {

        const int signalSize = aSignal->TimeSize();
        unsigned sampleIndex = 0;

        //Receiver Properties
        fDeltaT = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        fphiLO += 2. * LMCConst::Pi() * fLO_Frequency * fDeltaT;
        double tThisEventNSamples = fInterface->fTOld / fDeltaT;

    	std::vector<double> tKassParticleXP = fInterface->fTransmitter->ExtractParticleXP(fInterface->fTOld, true, fDeltaT, fInterface->fbWaveguide);
        double unitConversion = 1.;
        double excitationAmplitude = 0.;
        std::vector<double> tEFieldAtProbe;
        std::vector<double> dopplerFrequency;

        for (int mu=0; mu<fModeSet.size(); mu++)
		{
		    bool bTE = fModeSet[mu][0];
		    int l = fModeSet[mu][1];
		    int m = fModeSet[mu][2];
		    int n = fModeSet[mu][3];

		    std::vector<double> tE_normalized;
		    tE_normalized = fInterface->fField->GetNormalizedModeField(l,m,n,tKassParticleXP,1,bTE);
		    double cavityFIRSample = fFieldCalculator->GetCavityFIRSample(bTE, l, m, n, tKassParticleXP, fBypassTF).first;
		    dopplerFrequency = fInterface->fField->GetDopplerFrequency(l, m, n, tKassParticleXP);

		    double tAvgDotProductFactor = fInterface->fField->CalculateDotProductFactor(l, m, n, tKassParticleXP, tE_normalized, tThisEventNSamples);
		    double modeAmplitude = fInterface->fField->NormalizedEFieldMag(tE_normalized);

		    if (!fInterface->fbWaveguide)
		    {
		        // sqrt(4PIeps0) for Kass current si->cgs, sqrt(4PIeps0) for Jackson A_lambda coefficient cgs->si
		    	unitConversion = 1. / LMCConst::FourPiEps(); // see comment ^
		    	// Calculate propagating E-field with J \dot E.  cavityFIRSample units are [current]*[unitless].
		    	excitationAmplitude = tAvgDotProductFactor * modeAmplitude * cavityFIRSample * fInterface->fField->Z_TE(l,m,n,tKassParticleXP[7]) * 2. * LMCConst::Pi() / LMCConst::C() / 1.e2;
		    	tEFieldAtProbe = fInterface->fField->GetFieldAtProbe(l,m,n,1,tKassParticleXP,bTE);
		    }
		    else
		    {
		        // Waveguide default:  Use direct Kassiopeia power budget.
		    	if (fUseDirectKassPower)
		    	{
		    	    // replace signal amplitude with direct Kass power:
		    		unitConversion = 1.0;  // Kass power is already in Watts.
		    		std::vector<double> tTempKassParticleXP = {0.,0.,0.,0.,0.,0.,0.,tKassParticleXP[7],0.};
		    		double modeMax = fInterface->fField->GetNormalizedModeField(l,m,n,tTempKassParticleXP,0,1).back();
		    		double modeFrac = 0.; if (fabs(modeMax) > 0.) modeFrac = tE_normalized.back()/modeMax;
		    		excitationAmplitude = tAvgDotProductFactor*modeFrac*sqrt(tKassParticleXP[8]/2.);  // sqrt( modeFraction*LarmorPower/2 )
		    		tEFieldAtProbe = std::vector<double> {excitationAmplitude};
		    	}
		    	else
		    	{
		    	    // sqrt(4PIeps0) for Kass current si->cgs, sqrt(4PIeps0) for Jackson A_lambda coefficient cgs->si
		    		unitConversion = 1. / LMCConst::FourPiEps(); // see comment ^
		    		// Jackson |A| = 2piZ/c \int{J \dot E dV}, written assuming cgs units, and cavityFIRSample represents qvZ:
		    		excitationAmplitude = tAvgDotProductFactor * modeAmplitude * cavityFIRSample * 2. * LMCConst::Pi() / (LMCConst::C()*1.e2);
		    		// To extract the propagating E-field, scale the excitation amplitude as in the Pozar Poynting vector integral, p. 114:
		    		excitationAmplitude *= fInterface->fField->ScaleEPoyntingVector(tKassParticleXP[7]);
		    		tEFieldAtProbe = std::vector<double> {excitationAmplitude};
		    	}
		    } // Finished mode set.

		    for(int channelIndex = 0; channelIndex < fNChannels; ++channelIndex) // one channel per probe
		    {
		        sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample
		        // This scaling factor includes a 50 ohm impedance that is applied in signal processing, as well
		        // as other factors as defined above, e.g. 1/4PiEps0 if converting to/from c.g.s amplitudes.
		        double totalScalingFactor = sqrt(50.) * unitConversion;
		        fPowerCombiner->AddOneModeToCavityProbe(l, m, n, aSignal, tKassParticleXP, excitationAmplitude, tEFieldAtProbe[channelIndex], dopplerFrequency, fDeltaT, fphiLO, totalScalingFactor, sampleIndex, channelIndex, !(fInterface->fTOld > 0.) );
		    }
		}

        fInterface->fTOld += fDeltaT;
    	if (!fAliasingIsChecked)
    	{
    		if (!fOverrideAliasing)
    		{
    		    if (!CrossCheckAliasing( aSignal, dopplerFrequency[0] / 2. / LMCConst::Pi() ))
    		    {
    			    LERROR(lmclog,"There is an aliased HF frequency in the window.\n");
    			    return false;
    		    }
    		}
    		fAliasingIsChecked = true;
    	}

    	return true;
    }


    void CavitySignalGenerator::KassiopeiaInit(const std::string &aFile)
    {
        RunKassiopeia tRunKassiopeia;
        tRunKassiopeia.Run(aFile, fInterface);
        return;
    }


    void CavitySignalGenerator::WakeBeforeEvent()
    {
        fInterface->fPreEventCondition.notify_one();
        return;
    }

    bool CavitySignalGenerator::TryWakeAgain()
    {
    	int count = 0;
    	while (count < 10)
    	{
    		LPROG(lmclog,"Kass thread is unresponsive.  Trying again.\n");
    		LPROG( lmclog, "LMC about to try WakeBeforeEvent() again" );
    		std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    		WakeBeforeEvent();  // trigger Kass event.
    		if (!fInterface->fKassEventReady)  // Kass confirms event is underway.
    		{
    			return true;
    		}
    		count += 1;
    	}
 		return false;
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
        return (this->*fDoGenerateFunc)( aSignal );
    }


    bool CavitySignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

    bool CavitySignalGenerator::DoGenerateTime( Signal* aSignal )
    {
        ConfigureInterface( aSignal );
        RecordRunParameters( aSignal );

        if (fRandomPreEventSamples) RandomizeStartDelay();

        fPowerCombiner->SizeNChannels(fNChannels);
 	    if (fNChannels > 2)
 	    {
    		LERROR(lmclog,"The cavity simulation only supports up to 2 channels right now.");
        	throw std::runtime_error("Only 1 or 2 channels is allowed.");
        	return false;
 	    }

        int PreEventCounter = 0;

        if (fInterface->fTransmitter->IsKassiopeia())
        {
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
                    if (PreEventCounter > fNPreEventSamples)// finished pre-samples.  Start event.
                    {
                        fInterface->fPreEventInProgress = false;  // reset.
                        fInterface->fEventInProgress = true;
                        LPROG( lmclog, "LMC about to WakeBeforeEvent()" );
                        WakeBeforeEvent();  // trigger Kass event.
                    }
                }

                if (fInterface->fEventInProgress)  // fEventInProgress
                {

                    std::unique_lock< std::mutex >tLock( fInterface->fMutexDigitizer, std::defer_lock );



                    if (!fInterface->fKassEventReady)  // Kass confirms event is underway.
                    {

                    	fInterface->fSampleIndex = index; // 2-way trigger confirmation for Kass.
                    	tLock.lock();

                    	fInterface->fDigitizerCondition.wait( tLock );


                        if (fInterface->fEventInProgress)
                        {
                    		if (DriveMode(aSignal, index))
                    		{
                    			PreEventCounter = 0; // reset
                    		}
                    		else
                    		{
                    			fAliasedFrequencies = true;
                    			tLock.unlock();
                    			break;
                    		}
                        }
                        fInterface->fSampleIndex = index; // 2-way trigger confirmation for Kass.

                        tLock.unlock();


                	}
                 	else  // diagnose Kass
                 	{
                 		tLock.lock();
                 		std::this_thread::sleep_for(std::chrono::milliseconds(fThreadCheckTime));
                 		if (!fInterface->fKassEventReady)   // Kass event did start.  Continue but skip this sample.
                 		{
                 			tLock.unlock();
                 		}
                 		else    // Kass event has not started.
                 		{
                 			if ( fInterface->fEventInProgress )
                 			{
                 				if ( index < fNPreEventSamples+1 )  // Kass never started.
                 				{
                 					if (TryWakeAgain())
                 					{
                 						tLock.unlock();
                 					}
                 					else
                 					{
                     					LPROG(lmclog,"Locust is stopping because Kass has either stopped reporting, or never started.\n");
                     					tLock.unlock();
                 						break;
                 					}
                 				}
                 				else
                 				{
                 					LPROG(lmclog,"Locust infers that Kass has completed all events.\n");
                 					break;
                 				}
                 			}
                 			else
                 			{
                 				LWARN(lmclog, "Kass event terminated quickly.\n");
                 				tLock.unlock();
                 			}
                 		}
                 	} // diagnose Kass


                } // if fEventInProgress

            }  // for loop

            fInterface->fDoneWithSignalGeneration = true;
            LPROG( lmclog, "Finished signal loop." );
            fInterface->fWaitBeforeEvent = false;
            WakeBeforeEvent();
            tKassiopeia.join();  // finish thread

            if (fKassNeverStarted == true)
            {
            	throw std::runtime_error("Kassiopeia did not start.");
            	return false;
            }
            if (fAliasedFrequencies == true)
            {
            	throw std::runtime_error("There appears to be problematic HF aliasing in the window.  "
            			"See output above regarding \"Aliased frequency [] is below Nyquist frequency.\"  "
            			"Please try the unit test \"bin/testAliasingHF [-h]\" to optimize tuning.  Or, to override "
            			"this error please use this command line flag \"cavity-signal.override-aliasing\"=true ");
            	return false;
            }

        }  // fInterface->fTransmitter->IsKassiopeia()



    	return true;
    }



} /* namespace locust */


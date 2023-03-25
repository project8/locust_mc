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
        fNModes( 2 ),
        gxml_filename("blank.xml"),
        fphiLO(0.),
        fAvgDotProductFactor(0.),
        fNPreEventSamples( 150000 ),
        fThreadCheckTime(100),
        fKassNeverStarted( false ),
        fAliasedFrequencies( false ),
        fOverrideAliasing( false ),
        fOverrideStepsize( false ),
        fBypassTF( false ),
        fNormCheck( false ),
        fTE( true ),
        fIntermediateFile( false ),
        fUseDirectKassPower( false ),
        fAliasingIsChecked( false ),
        fInterface( new KassLocustInterface() )
    {
        KLInterfaceBootstrapper::get_instance()->SetInterface( fInterface );
    }

    CavitySignalGenerator::~CavitySignalGenerator()
    {
    }



    bool CavitySignalGenerator::ModeSelect(int l, int m, int n, bool eGun)
    {
    	if (eGun)
    	{
    		if (!fNormCheck)
    		{
    			if ((l==0)&&(m==1)&&(n==0))
    				return true;
    			else
    				return false;
    		}
    		else
    		{
    			if ((l<=fNModes)&&(m<=fNModes)&&(n<=fNModes))
    				return true;
    			else
    				return false;
    		}
    	}
    	else
    	{
    		if (!fNormCheck)
    		{
    			if ((l==0)&&(m==1)&&(n==1))
    				return true;
    			else
    				return false;
    		}
    		else
    		{
    			if ((l<=fNModes)&&(m<=fNModes)&&(n<=fNModes))
    				return true;
    			else
    				return false;
    		}
    	}
    }




    bool CavitySignalGenerator::Configure( const scarab::param_node& aParam )
    {

    	fTFReceiverHandler = new TFReceiverHandler;
    	if(!fTFReceiverHandler->Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver FIRHandler class");
    		exit(-1);
    		return false;
    	}

        if( aParam.has( "tf-receiver-filename" ) )
        {
            if (!fTFReceiverHandler->ReadHFSSFile())  // Read external file
            {
            	LERROR(lmclog,"FIR has not been generated.");
            	exit(-1);
            	return false;
            }
        }
        else // Generate analytic response function
        {
        	if ((aParam.has( "equivalent-circuit" ) ) && (aParam["equivalent-circuit"]().as_bool()))
        	{
        		fAnalyticResponseFunction = new EquivalentCircuit();
        		if ( !fAnalyticResponseFunction->Configure(aParam) )
        		{
        			LWARN(lmclog,"EquivalentCircuit was not configured.");
        			exit(-1);
        			return false;
        		}
        		else
        		{
        			if (!fTFReceiverHandler->ConvertAnalyticTFtoFIR(fAnalyticResponseFunction->GetInitialFreq(),fAnalyticResponseFunction->GetTFarray()))
        			{
        				LWARN(lmclog,"TF->FIR was not generated correctly.");
        				exit(-1);
        				return false;
        			}
        		}
        	}
        	else
        	{
        		fAnalyticResponseFunction = new DampedHarmonicOscillator();
        		if ( !fAnalyticResponseFunction->Configure(aParam) ||
        				(!CrossCheckCavityConfig()) )
        		{
        			LERROR(lmclog,"DampedHarmonicOscillator was not configured.");
        			exit(-1);
        			return false;
        		}
        		if (!fTFReceiverHandler->ConvertAnalyticGFtoFIR(fAnalyticResponseFunction->GetGFarray()))
        		{
        			LERROR(lmclog,"GF->FIR was not generated.");
        			exit(-1);
        			return false;
        		}
        	}
        } // aParam.has( "tf-receiver-filename" )
        fdtFilter = fTFReceiverHandler->GetFilterResolution();


        if( aParam.has( "e-gun" ) )
        {
        	fInterface->fE_Gun = aParam["e-gun"]().as_bool();
        }

        if (fInterface->fE_Gun)
        {
        	fInterface->fField = new RectangularWaveguide;
        	fPowerCombiner = new WaveguideModes;
        	if ( aParam.has( "waveguide-short" ) )
        	{
        		fPowerCombiner->SetWaveguideShortIsPresent(aParam["waveguide-short"]().as_bool());
        	}
        	else
        	{
        		// This is the same as the default case in LMCPowerCombiner
        		fPowerCombiner->SetWaveguideShortIsPresent( true );
        	}
        	// Allow back reaction only if short is present:
        	if ( fPowerCombiner->GetWaveguideShortIsPresent() )
        	{
        		if ( aParam.has( "back-reaction" ) )
        		{
        			fInterface->fBackReaction = aParam["back-reaction"]().as_bool();
        		}
        		else
        		{
        			fInterface->fBackReaction = true;
        		}
        	}
        }
        else
        {
        	fInterface->fField = new CylindricalCavity;
        	fPowerCombiner = new CavityModes;
        	if ( aParam.has( "back-reaction" ) )
        	{
        		fInterface->fBackReaction = aParam["back-reaction"]().as_bool();
                if( aParam.has( "override-stepsize" ) )
                {
                    fOverrideStepsize = aParam["override-stepsize"]().as_bool();
                }
        		if (( !fInterface->fBackReaction ) && (fOverrideStepsize == false))
        		{
        			LERROR(lmclog,"If \"back-reaction\"=false, please check that the Kass "
        					"stepsize is 0.135 orbits or smaller.  Otherwise there can be "
        					"skipped digitizer samples when running in a multi-core cluster "
        					"environment.  This is a to-do item to be fixed with a planned 2-step trigger. "
        					"To override this message, use this command-line flag:  \"cavity-signal.override-stepsize\"=true");
        			exit(-1);
        			return false;
        		}
        	}
        }
        if( aParam.has( "n-modes" ) )
        {
        	fNModes = aParam["n-modes"]().as_int();
        }

        fPowerCombiner->SetNCavityModes(fNModes);
        if(!fPowerCombiner->Configure(aParam))
		{
			LERROR(lmclog,"Error configuring PowerCombiner.");
			exit(-1);
			return false;
		}


        fInterface->fField->SetNModes(fNModes);
        if (!fInterface->fField->Configure(aParam))
        {
        	LERROR(lmclog,"Error configuring LMCField.");
        	exit(-1);
        	return false;
        }


        if( aParam.has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam["lo-frequency"]().as_double();
        }

        if( aParam.has( "event-spacing-samples" ) )
        {
            fNPreEventSamples = aParam["event-spacing-samples"]().as_int();
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
        if( aParam.has( "te-modes" ) )
        {
        	fTE = aParam["te-modes"]().as_bool();
        }
        if( aParam.has( "intermediate-file" ) )
        {
        	fIntermediateFile = aParam["intermediate-file"]().as_bool();
        }
        if( aParam.has( "direct-kass-power" ) )
        {
        	fUseDirectKassPower = aParam["direct-kass-power"]().as_bool();
        }
        if( aParam.has( "xml-filename" ) )
        {
            gxml_filename = aParam["xml-filename"]().as_string();
        }

        if( aParam.has( "transmitter" ))
        {
        	int ntransmitters = 0;

        	if(aParam["transmitter"]().as_string() == "kass-current")
        	{
        		ntransmitters += 1;
        		fInterface->fTransmitter = new KassCurrentTransmitter;
        		if(!fInterface->fTransmitter->Configure(aParam))
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


        fFieldCalculator = new FieldCalculator();
        if(!fFieldCalculator->Configure(aParam))
        {
            LERROR(lmclog,"Error configuring receiver FieldCalculator class from CavitySignal.");
        }
        fInterface->fConfigureKass = new ConfigureKass();
        fInterface->fConfigureKass->SetParameters( aParam );


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

    	CavityUtility aCavityUtility;
    	double timeResolution = fAnalyticResponseFunction->GetDHOTimeResolution();
    	double thresholdFactor = fAnalyticResponseFunction->GetDHOThresholdFactor();
    	double cavityFrequency = fAnalyticResponseFunction->GetCavityFrequency();
    	double qExpected = fAnalyticResponseFunction->GetCavityQ();
    	if (!aCavityUtility.CheckCavityQ( timeResolution, thresholdFactor, cavityFrequency, qExpected ))
    	{
        	LERROR(lmclog,"The cavity Q does not look quite right.  Please tune the configuration "
        			"with the unit test as in bin/testLMCCavity [-h]");
    		return false;
    	}
    	else
    	{
    		return true;
    	}
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






    double CavitySignalGenerator::ScaleEPoyntingVector(double fcyc)
    {
    	// This function calculates the coefficients of the Poynting vector integral
    	// in the TE10 mode in WR42.  It then returns the sqrt of the half of the propagating
    	// power that is moving toward the antenna.
    	// After Pozar p. 114:
    	double k = fcyc / LMCConst::C();
    	double k1 = LMCConst::Pi() / fInterface->fField->GetDimX();
    	double beta = sqrt( k*k - k1*k1 );
    	double areaIntegral = fcyc * LMCConst::MuNull() * pow(fInterface->fField->GetDimX(),3.) * fInterface->fField->GetDimY() * beta / 4. / LMCConst::Pi() / LMCConst::Pi();
    	// sqrt of propagating power gives amplitude of E
    	return sqrt(areaIntegral);
    }

    bool CavitySignalGenerator::DriveMode(Signal* aSignal, unsigned index)
    {

        const int signalSize = aSignal->TimeSize();
        unsigned sampleIndex = 0;
        const unsigned nChannels = fNChannels;

        //Receiver Properties
        fDeltaT = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        fphiLO += 2. * LMCConst::Pi() * fLO_Frequency * fDeltaT;
        double tThisEventNSamples = fInterface->fTOld / fDeltaT;

    	std::vector<double> tKassParticleXP = fInterface->fTransmitter->ExtractParticleXP(fInterface->fTOld, true, fDeltaT, fInterface->fE_Gun);
        double unitConversion = 1.;
        double excitationAmplitude = 0.;
        double tEFieldAtProbe = 0.;
        std::vector<double> dopplerFrequency;


    	for (int l=0; l<fNModes; l++)
    	{
    		for (int m=1; m<fNModes; m++)
    		{
    			for (int n=0; n<fNModes; n++)
    			{
    				if (ModeSelect(l, m, n, fInterface->fE_Gun))
    				{
    					std::vector<double> tE_normalized;
    					tE_normalized = fInterface->fField->GetNormalizedModeField(l,m,n,tKassParticleXP);
    					double cavityFIRSample = fFieldCalculator->GetCavityFIRSample(tKassParticleXP, fBypassTF).first;
    					dopplerFrequency = fInterface->fField->GetDopplerFrequency(l, m, n, tKassParticleXP);
    					fAvgDotProductFactor = 1. / ( tThisEventNSamples + 1 ) * ( fAvgDotProductFactor * tThisEventNSamples + fInterface->fField->GetDotProductFactor(tKassParticleXP, tE_normalized, fIntermediateFile) );  // unit velocity \dot unit theta
    					double modeAmplitude = tE_normalized.back();


    					if (!fInterface->fE_Gun)
    					{
    						// sqrt(4PIeps0) for Kass current si->cgs, sqrt(4PIeps0) for Jackson A_lambda coefficient cgs->si
    						unitConversion = 1. / LMCConst::FourPiEps(); // see comment ^
    						// Calculate propagating E-field with J \dot E.  cavityFIRSample units are [current]*[unitless].
    						excitationAmplitude = fAvgDotProductFactor * modeAmplitude * cavityFIRSample * fInterface->fField->Z_TE(l,m,n,tKassParticleXP[7]) * 2. * LMCConst::Pi() / LMCConst::C() / 1.e2;
    						tEFieldAtProbe = fInterface->fField->GetFieldAtProbe(l,m,n,1);
    					}
    					else
    					{
    						// sqrt(4PIeps0) for Kass current si->cgs, sqrt(4PIeps0) for Jackson A_lambda coefficient cgs->si
    						unitConversion = 1. / LMCConst::FourPiEps(); // see comment ^
    						// Calculate propagating E-field with J \dot E and integrated Poynting vector.  cavityFIRSample units are [current]*[ohms].
    						excitationAmplitude = fAvgDotProductFactor * modeAmplitude * ScaleEPoyntingVector(tKassParticleXP[7]) * cavityFIRSample * 2. * LMCConst::Pi() / LMCConst::C() / 1.e2;

    						// Optional cross-check:  Use direct Kassiopeia power budget.  Assume x_electron = 0.
    						if (fUseDirectKassPower)
    						{
    							// override one-way signal amplitude with direct Kass power:
    							unitConversion = 1.0;  // Kass power is already in Watts.
        						excitationAmplitude = fAvgDotProductFactor*sqrt(tKassParticleXP[8]/2.);  // sqrt( modeFraction*LarmorPower/2 )
    						}

    					}

    					for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)  // one channel per probe.
    					{
    						sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample

    						// This scaling factor includes a 50 ohm impedance that applied in signal processing, as well
    						// as other factors as defined above, e.g. 1/4PiEps0 if converting to/from c.g.s amplitudes.
    						double totalScalingFactor = sqrt(50.) * unitConversion;
    						fPowerCombiner->AddOneModeToCavityProbe(aSignal, tKassParticleXP, excitationAmplitude, tEFieldAtProbe, dopplerFrequency, fDeltaT, fphiLO, totalScalingFactor, sampleIndex, !(fInterface->fTOld > 0.) );
    						if (fNormCheck) fPowerCombiner->AddOneSampleToRollingAvg(l, m, n, excitationAmplitude, sampleIndex);
    					}

    				} // ModeSelect
    			} // n
    		} // m
    	} // l

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

        int PreEventCounter = 0;
        fFieldCalculator->SetNFilterBinsRequired( 1. / (fAcquisitionRate*1.e6*aSignal->DecimationFactor()) );
        fFieldCalculator->SetFilterSize( fTFReceiverHandler->GetFilterSize() );

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
                 				if ( index < fNPreEventSamples+1 )  // Kass never started at all.
                 				{
                 					LERROR(lmclog,"Kass thread is unresponsive.  Exiting.\n");
                 					fKassNeverStarted = true;
                 				}
                 				tLock.unlock();   // Kass either started or not, but is now finished.
                 				break;
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


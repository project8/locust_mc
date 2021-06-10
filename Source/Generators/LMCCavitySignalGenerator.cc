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
        fLO_Frequency( 0.),
		fNModes( 1 ),
        gxml_filename("blank.xml"),
        fphiLO(0.),
		fNPreEventSamples( 150000 ),
		fThreadCheckTime(100),
		fKassNeverStarted( false ),
		fSkippedSamples( false ),
		fInterface( new KassLocustInterface() )
    {
        fRequiredSignalState = Signal::kFreq;
        KLInterfaceBootstrapper::get_instance()->SetInterface( fInterface );
    }

    CavitySignalGenerator::~CavitySignalGenerator()
    {
    }


    void CavitySignalGenerator::ReadBesselZeroes(std::string filename, bool prime)
    {
        std::ifstream input( filename );
        int n = 0; int k = 0; double zero = 0.;
        for( std::string line; getline( input, line ); )
        {
            std::stringstream ss(line);
            ss >> n;
            ss >> k;
            ss >> zero;
            if (prime)
            {
            	fInterface->fBesselNKPrimeZeros.resize(n+1);
            	fInterface->fBesselNKPrimeZeros[n].resize(k+1);
            	fInterface->fBesselNKPrimeZeros[n][k] = zero;
            }
            else
            {
            	fInterface->fBesselNKZeros.resize(n+1);
            	fInterface->fBesselNKZeros[n].resize(k+1);
            	fInterface->fBesselNKZeros[n][k] = zero;
            }

//            printf("zero is %g and zero01 is %g\n\n", fInterface->fBesselNKPrimeZeros[n][k], fInterface->fBesselNKPrimeZeros[0][1]);
        }

    }


    void CavitySignalGenerator::CheckNormalization()
    {
    	CylindricalCavity aCavity;

    	printf("\\epsilon\\int{|E_xlm|^2 dV} = \\mu\\int{|H_xlm|^2 dV} ?\n\n");
    	for (int l=0; l<3; l++)
    		for (int m=1; m<4; m++)
    			for (int n=1; n<4; n++)
    			{
    		    	printf("TE%d%d%d E %g H %g\n", l, m, n, LMCConst::EpsNull()*aCavity.Integrate(l,m,n,1,1), LMCConst::MuNull()*aCavity.Integrate(l,m,n,1,0));
    			}
    	for (int l=0; l<3; l++)
    		for (int m=1; m<3; m++)
    			for (int n=1; n<3; n++)
    			{
//    		    	printf("TM%d%d%d E %g H %g\n", l, m, n, LMCConst::EpsNull()*aCavity.Integrate(l,m,n,0,1), LMCConst::MuNull()*aCavity.Integrate(l,m,n,0,0));
    			}
    }



    void CavitySignalGenerator::PrintModeMaps()
    {
    	char buffer[60];
    	CylindricalCavity aCavity;

    	for (int l=0; l<5; l++)
    		for (int m=1; m<5; m++)
    			for (int n=1; n<5; n++)
    			{
    				printf("l m n is %d %d %d\n", l, m, n);
    				double tNormalizationTE_E = pow(aCavity.Integrate(l,m,n,1,1),0.5);
    				int a = sprintf(buffer, "output/ModeMapTE%d%d%d_E.txt", l, m, n);
    				const char *fpname = buffer;
    				FILE *fpTE_E = fopen(fpname, "w");
    				for (unsigned i=0; i<fInterface->fnPixels; i++)
    				{
    					double r = (double)i/fInterface->fnPixels*fInterface->fR;
    					for (unsigned j=0; j<fInterface->fnPixels; j++)
    					{
    						double theta = (double)j/fInterface->fnPixels*2.*LMCConst::Pi();
    						std::vector<double> tTE_E = aCavity.TE_E(l,m,n,r,theta,0.05);
    						fprintf(fpTE_E, "%10.4g %10.4g %10.4g %10.4g\n", r, theta, tTE_E.front()/tNormalizationTE_E, tTE_E.back()/tNormalizationTE_E);
    					}
    				}
    				fclose (fpTE_E);
    			}
    }



    bool CavitySignalGenerator::Configure( const scarab::param_node& aParam )
    {

    	if(!fInterface->fTFReceiverHandler.Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver FIRHandler class");
    	}

        if(!fInterface->fTFReceiverHandler.ReadHFSSFile())
        {
            return false;
        }

        if( aParam.has( "cavity-radius" ) )
        {
            fInterface->fR = aParam["cavity-radius"]().as_double();
            fInterface->fR = 500.;
        }

        if( aParam.has( "cavity-length" ) )
        {
            fInterface->fL = aParam["cavity-length"]().as_double();
        }


        fInterface->dtFilter = fInterface->fTFReceiverHandler.GetFilterResolution();
    	FieldBuffer aFieldBuffer;
    	fInterface->eCurrentBuffer = aFieldBuffer.InitializeBuffer(1, 1, fInterface->fTFReceiverHandler.GetFilterSize());


        scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
        ReadBesselZeroes((dataDir / "BesselZeros.txt").string(), 0 );
        ReadBesselZeroes((dataDir / "BesselPrimeZeros.txt").string(), 1 );
        CheckNormalization();
        PrintModeMaps();

        if( aParam.has( "transmitter" ))
        {
        	int ntransmitters = 0;

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
        		LERROR(lmclog,"LMCCavitySignalGenerator needs a single transmitter.  Please choose transmitter:kassiopeia in the config file.");
                exit(-1);
        	}
        }
        else
        {
    		LERROR(lmclog,"LMCCavitySignalGenerator has been configured without a transmitter.  Please choose transmitter:antenna or transmitter:planewave or transmitter:kassiopeia in the config file.");
            exit(-1);
        }

        if( aParam.has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam["lo-frequency"]().as_double();
        }

        if( aParam.has( "n-modes" ) )
        {
            fNModes = aParam["n-modes"]().as_int();
        }

        if( aParam.has( "event-spacing-samples" ) )
        {
            fNPreEventSamples = aParam["event-spacing-samples"]().as_int();
        }
        if( aParam.has( "thread-check-time" ) )
        {
            fThreadCheckTime = aParam["thread-check-time"]().as_int();
        }
        if( aParam.has( "xml-filename" ) )
        {
            gxml_filename = aParam["xml-filename"]().as_string();
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


    bool CavitySignalGenerator::DoGenerateTimeKass( Signal* aSignal )
    {


        //n samples for event spacing in Kass.
        int PreEventCounter = 0;
        fInterface->nFilterBinsRequired = (int) std::min( 1. / ((fAcquisitionRate*1.e6*aSignal->DecimationFactor()) / fInterface->dtFilter), (double) fInterface->fTFReceiverHandler.GetFilterSize());


        if (!fTransmitter->IsKassiopeia())
        {
        	for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        	{
//        		DriveAntenna(fp, PreEventCounter, index, aSignal, nfilterbins, dtfilter);
        	}  // for loop
        	return true;
        }


        if (fTransmitter->IsKassiopeia())
        {
        	bool fTruth = false;
        	int startingIndex;
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
/*
                      		if (DriveAntenna(fp, startingIndex, index, aSignal, nfilterbins, dtfilter))
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
*/
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

    bool CavitySignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
    	const unsigned nchannels = 1;
    	double fLO_frequency = 20.03e9;
    	double fRF_frequency = 20.0e9;
    	double fAmplitude = 1.e-9;
        double LO_phase = 0.;
        double voltage_phase = 0.;

        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
            for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
            {
                LO_phase = 2.*LMCConst::Pi()*fLO_frequency*(double)index/(fAcquisitionRate*1.e6);
                voltage_phase = 2.*LMCConst::Pi()*fRF_frequency*(double)index/(fAcquisitionRate*1.e6);

                if (index == aSignal->FreqSize()/6 ) // pick a frequency bin and put some signal in it.
                {
                	aSignal->SignalFreqComplex()[ch*aSignal->TimeSize() + index][0] += sqrt(50.)*fAmplitude;
                	aSignal->SignalFreqComplex()[ch*aSignal->TimeSize() + index][1] += sqrt(50.)*fAmplitude;
                }
            }
        }

        aSignal->ToState(Signal::kTime);

        return true;
    }

    bool CavitySignalGenerator::DoGenerateTime( Signal* aSignal )
    {

    	const unsigned nchannels = 1;
    	double fLO_frequency = 20.e9;
    	double fRF_frequency = 20.05e9;
    	double fAmplitude = 1.e-8;
        double LO_phase = 0.;
        double voltage_phase = 0.;

        for (unsigned ch = 0; ch < nchannels; ++ch)
        {
            for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
            {
                LO_phase = 2.*LMCConst::Pi()*fLO_frequency*(double)index/(fAcquisitionRate*1.e6);
                voltage_phase = 2.*LMCConst::Pi()*fRF_frequency*(double)index/(fAcquisitionRate*1.e6);

                aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][0] += sqrt(50.)*fAmplitude*cos(voltage_phase-LO_phase);
                aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][1] += sqrt(50.)*fAmplitude*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);
            }
        }

        return true;
    }




} /* namespace locust */


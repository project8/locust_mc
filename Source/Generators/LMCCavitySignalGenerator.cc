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
		fThreadCheckTime(200),
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

    std::pair<double, double> CavitySignalGenerator::TE(int l, int m, int n, double r, double theta, double z)
    {
    	// from "Techniques of Microwave Measurements", C. G. Montgomery
    	double x_lm = fBesselNKPrimeZeros[l][m];
    	double k1 = 2.*x_lm / (2.*fR);
    	double k3 = n * LMCConst::Pi() / fL; // n = 1
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double tEr = -l * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
    	double tEtheta = -jPrime * cos(l*theta) * sin(k3*z);
        return std::make_pair(tEr, tEtheta);
    }


    void CavitySignalGenerator::ReadFile(std::string filename, std::vector<std::vector<double> > &data)
    {
        std::ifstream input( filename );
        int n = 0; int k = 0; double zero = 0.;
        for( std::string line; getline( input, line ); )
        {
            std::stringstream ss(line);
            ss >> n;
            ss >> k;
            ss >> zero;
            data.resize(n+1);
            data[n].resize(k+1);
            data[n][k] = zero;
//            printf("zero is %g and zero01 is %g\n\n", data[n][k], data[0][1]);
        }

    }

    void CavitySignalGenerator::PrintModeMap()
    {
    	FILE *fp = fopen("ModeMapDebug.txt", "w");
    	fR = 0.1; //m, parametrize me.
    	fL = 0.2; // m, parametrize me.
    	fnPixels = 100; // me too.
    	for (unsigned tR=0; tR<fnPixels; tR++)
    	{
    		double r = fR/fnPixels*tR;
    		for (unsigned tTheta=0; tTheta<fnPixels; tTheta++)
    		{
    			double theta = 2.*LMCConst::Pi()/fnPixels*tTheta;
    			std::pair<double, double> tFieldValues = TE(0,1,1,r,theta,0.1);
             	fprintf(fp, "%10.4g %10.4g %10.4g %10.4g\n", r, theta, tFieldValues.first, tFieldValues.second);
    		}
    	}
    	fclose (fp);


    }



    bool CavitySignalGenerator::Configure( const scarab::param_node& aParam )
    {
        scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
        ReadFile((dataDir / "BesselZeros.txt").string(), fBesselNKZeros );
        ReadFile((dataDir / "BesselPrimeZeros.txt").string(), fBesselNKPrimeZeros );
        PrintModeMap();

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

/*    void CavitySignalGenerator::InitializeFieldPoints(std::vector< Channel<Receiver*> > allRxChannels)
    {
	for(int channelIndex = 0; channelIndex < fNChannels; ++channelIndex)
	{
            for(int elementIndex = 0; elementIndex < fNElementsPerStrip; ++elementIndex)
            {
            	fTransmitter->InitializeFieldPoint(allRxChannels[channelIndex][elementIndex]->GetPosition());
            }
	}

    }
    */

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

//        InitializeFieldPoints(allRxChannels);

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


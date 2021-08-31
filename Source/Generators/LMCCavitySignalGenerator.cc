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
		fNModes( 3 ),
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

    std::vector<int> CavitySignalGenerator::ModeFilter(unsigned whichMode)
    {
    	// This is not really a mode filter, but could be modifed/replaced to act like one.

    	unsigned modeCounter = 0;
    	std::vector<int> modeIndex;
    	modeIndex.resize(3);

    	for (int l=0; l<fNModes; l++)
    	{
    		modeIndex[0] = l;
    		for (int m=1; m<fNModes+1; m++)
    		{
        		modeIndex[1] = m;
    			for (int n=1; n<fNModes+1; n++)
    			{
            		modeIndex[2] = n;
    	    		modeCounter += 1;
//    	    		printf("l is %d, m is %d, n is %d, modeCounter is %d\n", l, m, n, modeCounter);
    	    		if (modeCounter > whichMode)
    	    		{
    	    	    	return modeIndex;
    	    		}
    			}
    		}
    	}


    }


    std::vector<std::vector<std::vector<double>>> CavitySignalGenerator::CalculateNormFactors(int nModes, bool TE)
    {

        LPROG( lmclog, "Calculating mode normalization factors ... " );

    	std::vector<std::vector<std::vector<double>>> aModeNormFactor;
    	aModeNormFactor.resize(nModes);

    	for (unsigned m=0; m<nModes; m++)
    	{
    		aModeNormFactor[m].resize(nModes);
        	for (unsigned n=0; n<nModes; n++)
        	{
        		aModeNormFactor[m][n].resize(nModes);
        	}
    	}


    	for (unsigned l=0; l<nModes; l++)
    	{
        	for (unsigned m=1; m<nModes; m++)
        	{
            	for (unsigned n=1; n<nModes; n++)
            	{

            		if (TE)
            		{
            			aModeNormFactor[l][m][n] = 1./fInterface->fField->Integrate(l,m,n,1,1)/LMCConst::EpsNull();
            		}
            		else
            		{
            			aModeNormFactor[l][m][n] = 1./fInterface->fField->Integrate(l,m,n,0,1)/LMCConst::EpsNull();
            		}
            	}
        	}
    	}

    	return aModeNormFactor;
    }


    void CavitySignalGenerator::CheckNormalization()
    {

    	printf("\n \\epsilon\\int{|E_xlm|^2 dV} = \\mu\\int{|H_xlm|^2 dV} ?\n\n");
    	for (int l=0; l<fNModes; l++)
    	{
    		for (int m=1; m<fNModes; m++)
    		{
    			for (int n=1; n<fNModes; n++)
    			{
    				double normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];
       		    	printf("TE%d%d%d E %.4f H %.4f\n", l, m, n, LMCConst::EpsNull()*fInterface->fField->Integrate(l,m,n,1,1)*normFactor,
        		    		LMCConst::MuNull()*fInterface->fField->Integrate(l,m,n,1,0)*normFactor);
    			}
    		}
    	}


    	for (int l=0; l<fNModes; l++)
    	{
    		for (int m=1; m<fNModes; m++)
    		{
    			for (int n=1; n<fNModes; n++)
    			{
    				double normFactor = fInterface->fField->GetNormFactorsTM()[l][m][n];
    		    	printf("TM%d%d%d E %.4f H %.4f\n", l, m, n, LMCConst::EpsNull()*fInterface->fField->Integrate(l,m,n,0,1)*normFactor,
    		    			LMCConst::MuNull()*fInterface->fField->Integrate(l,m,n,0,0)*normFactor);
    			}
    		}
    	}

    	printf("The modes normalized as above are available for use in the simulation.\n");
    }



    void CavitySignalGenerator::PrintModeMaps()
    {
    	char buffer[60];
    	unsigned modeCounter = 0;

    	for (int l=0; l<5; l++)
    		for (int m=1; m<5; m++)
    			for (int n=1; n<5; n++)
    			{
    				printf("l m n is %d %d %d\n", l, m, n);
    				double tNormalizationTE_E = fInterface->fField->GetNormFactorsTE()[l][m][n];
    				int a = sprintf(buffer, "output/ModeMapTE%d%d%d_E.txt", l, m, n);
    				const char *fpname = buffer;
    				FILE *fpTE_E = fopen(fpname, "w");
    				for (unsigned i=0; i<fInterface->fnPixels; i++)
    				{
    					double r = (double)i/fInterface->fnPixels*fInterface->fR;
    					for (unsigned j=0; j<fInterface->fnPixels; j++)
    					{
    						double theta = (double)j/fInterface->fnPixels*2.*LMCConst::Pi();
    						std::vector<double> tTE_E = fInterface->fField->TE_E(l,m,n,r,theta,0.0);
    						fprintf(fpTE_E, "%10.4g %10.4g %10.4g %10.4g\n", r, theta, tTE_E.front()*tNormalizationTE_E, tTE_E.back()*tNormalizationTE_E);
    					}
    				}
    				fclose (fpTE_E);
    				modeCounter += 1;
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
        }

        if( aParam.has( "cavity-length" ) )
        {
            fInterface->fL = aParam["cavity-length"]().as_double();
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


		fPowerCombiner = new CavityModes;
		if(!fPowerCombiner->Configure(aParam))
		{
			LERROR(lmclog,"Error configuring CavityModes:PowerCombiner.");
			exit(-1);
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
    		LERROR(lmclog,"LMCCavitySignalGenerator has been configured without a transmitter.  Please choose transmitter:antenna or transmitter:planewave or transmitter:kassiopeia in the config file.");
            exit(-1);
        }

        fInterface->dtFilter = fInterface->fTFReceiverHandler.GetFilterResolution();
    	FieldBuffer aFieldBuffer;
    	fInterface->eCurrentBuffer = aFieldBuffer.InitializeBuffer(1, 1, fInterface->fTFReceiverHandler.GetFilterSize());
    	fInterface->fField = new CylindricalCavity;


        scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
        ReadBesselZeroes((dataDir / "BesselZeros.txt").string(), 0 );
        ReadBesselZeroes((dataDir / "BesselPrimeZeros.txt").string(), 1 );


        fInterface->fField->SetNormFactorsTE(CalculateNormFactors(fNModes,1));
        fInterface->fField->SetNormFactorsTM(CalculateNormFactors(fNModes,0));
        CheckNormalization();  // E fields integrate to 1.0 for both TE and TM modes.  H fields near 1.0
//        PrintModeMaps();

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


    double CavitySignalGenerator::GetModeScalingFactor(std::vector<double> tKassParticleXP, int channelIndex)
    {
    	// Scale excitation amplitude at probe according to mode map:
    	std::vector<double> tTE_E_electron = fInterface->fField->TE_E(0,1,1,tKassParticleXP[0],0.,tKassParticleXP[2]);
    	double thetaProbe = tKassParticleXP[1] + fPowerCombiner->GetCavityProbeTheta()[channelIndex];
    	std::vector<double> tTE_E_probe = fInterface->fField->TE_E(0,1,1,fInterface->fR*0.95,thetaProbe,fPowerCombiner->GetCavityProbeZ()[channelIndex]);
        double modeAmplitudeFactor = tTE_E_probe.back() / tTE_E_electron.back();
        //        	printf("modeamplitudeFactor is %g / %g = %g\n", tTE_E_probe.back(), tTE_E_electron.back(), modeAmplitudeFactor); getchar();

    	return modeAmplitudeFactor;
    }


    double CavitySignalGenerator::GetCavityDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized)
    {
    	double tThetaParticle = tKassParticleXP[1];
    	double tEtheta = aTE_E_normalized.back();
    	double tEx = -sin(tThetaParticle) * tEtheta;
    	double tEy = cos(tThetaParticle) * tEtheta;
    	double tEmag = fabs(tEtheta);
    	// fix-me:  should below be instantaneous velocity or center of motion?
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double tVmag = pow(tVx*tVx + tVy*tVy, 0.5);
//    	printf("tEmag is %g, r is %g, and theta is %g\n", tEmag, tKassParticleXP[0], tKassParticleXP[1]); getchar();
    	return fabs(tEx*tVx + tEy*tVy)/tEmag/tVmag;  // fabs ( unit J \cdot unit E )
    }

    std::vector<double> CavitySignalGenerator::GetCavityNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP)
     {
     	double tR = tKassParticleXP[0];
     	double tZ = tKassParticleXP[2];
     	std::vector<double> tTE_E_electron = fInterface->fField->TE_E(l,m,n,tR,0.,tZ);
 		double normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];  // select mode 0,1,1

 		auto it = tTE_E_electron.begin();
 		while (it != tTE_E_electron.end())
 		{
 			if (!isnan(*it))
 				(*it) *= normFactor;
 			*it++;
 		}
     	return tTE_E_electron;  // return normalized field.
     }

    double CavitySignalGenerator::GetCavityFIRSample(std::vector<double> tKassParticleXP, std::vector<std::deque<double>> aLocalFIRfrequencyBuffer, std::vector<std::deque<double>> aLocalElementFIRBuffer,int nFilterBinsRequired, double dtFilter)
    {
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double orbitPhase = tKassParticleXP[6];  // radians
    	double fieldFrequency = tKassParticleXP[7];  // rad/s
    	double vMag = pow(tVx*tVx + tVy*tVy,0.5);
    	double convolution = 0.0;

		// populate FIR filter with frequency for just this sample interval:
		for (int i=0; i < nFilterBinsRequired; i++)
		{
			aLocalFIRfrequencyBuffer[0].push_back(fieldFrequency);  // rad/s
			aLocalFIRfrequencyBuffer[0].pop_front();
		}
		std::deque<double>::iterator it = aLocalFIRfrequencyBuffer[0].begin();
		while (it != aLocalFIRfrequencyBuffer[0].end())
		{
			orbitPhase += (*it)*dtFilter;

			if (*it != 0.)
			{
				aLocalElementFIRBuffer[0].push_back(LMCConst::Q()*vMag*cos(orbitPhase));
			}
			else
			{
				aLocalElementFIRBuffer[0].push_back(0.);
			}
			aLocalElementFIRBuffer[0].pop_front();

			*it++;
		}

		convolution=fInterface->fTFReceiverHandler.ConvolveWithFIRFilter(aLocalElementFIRBuffer[0]);
		return convolution;

    }


    bool CavitySignalGenerator::DriveMode(Signal* aSignal, int nFilterBinsRequired, double dtFilter, unsigned index)
    {
        const int signalSize = aSignal->TimeSize();
        unsigned sampleIndex = 0;
        const unsigned nChannels = fNChannels;

        //Receiver Properties
        fphiLO += 2. * LMCConst::Pi() * fLO_Frequency * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());


    	std::vector<double> tKassParticleXP = fInterface->fTransmitter->ExtractParticleXP(fInterface->fTOld);
    	std::vector<double> tTE_E_normalized = GetCavityNormalizedModeField(0,1,1,tKassParticleXP); // lmn 011 for now.
        double dotProductFactor = GetCavityDotProductFactor(tKassParticleXP, tTE_E_normalized);  // unit velocity \dot unit theta
    	double modeAmplitude = tTE_E_normalized.back();  // absolute E_theta at electron

// fix-me:  We may need more precise nFilterBinsRequired.
    	std::vector<std::deque<double>> tLocalFIRfrequencyBuffer = fInterface->FIRfrequencyBufferCopy;  // copy from Kass buffer.
    	std::vector<std::deque<double>> tLocalElementFIRBuffer = fInterface->ElementFIRBufferCopy;
        double tFirSample = GetCavityFIRSample(tKassParticleXP, tLocalFIRfrequencyBuffer, tLocalElementFIRBuffer, fInterface->nFilterBinsRequired, fInterface->dtFilter);
        std::vector<std::deque<double>>().swap(tLocalFIRfrequencyBuffer);  // release memory
        std::vector<std::deque<double>>().swap(tLocalElementFIRBuffer);

        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)  // one channel per probe.
        {
        	sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample
//       	double modeScalingFactor = GetModeScalingFactor(tKassParticleXP, channelIndex);  // scale over to the probe?
        	double modeScalingFactor = 1.0; // override for now.
        	double totalScalingFactor = modeScalingFactor * dotProductFactor * modeAmplitude;
        	fPowerCombiner->AddOneModeToCavityProbe(aSignal, tFirSample, fphiLO, totalScalingFactor, fPowerCombiner->GetCavityProbeImpedance(), sampleIndex);
        }

        fInterface->fTOld += 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

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

        // dump to text file for debugging.
    	fp = fopen("currentValues.txt", "w");



        int PreEventCounter = 0;
		fPowerCombiner->SetCavityProbeLocations(fNChannels, fInterface->fL);
        int nFilterBins = fInterface->fTFReceiverHandler.GetFilterSize();
        InitializeBuffers(nFilterBins);


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
                    		if (DriveMode(aSignal, fInterface->nFilterBinsRequired, fInterface->dtFilter, index))
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

        }  // fInterface->fTransmitter->IsKassiopeia()

    	fclose (fp);


    	return true;
    }


    void CavitySignalGenerator::InitializeBuffers(unsigned filterbuffersize)
    {
    	FieldBuffer aFieldBuffer;
    	fInterface->ElementFIRBuffer = aFieldBuffer.InitializeBuffer(1, 1, filterbuffersize);
    	fInterface->ElementFIRBufferCopy = aFieldBuffer.InitializeBuffer(1, 1, filterbuffersize);
    	fInterface->FIRfrequencyBuffer = aFieldBuffer.InitializeBuffer(1, 1, filterbuffersize);
    	fInterface->FIRfrequencyBufferCopy = aFieldBuffer.InitializeBuffer(1, 1, filterbuffersize);
    }





} /* namespace locust */


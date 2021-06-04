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
    	fR(0.1), //m
    	fL(0.2), // m
    	fnPixels(100),
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

    double CavitySignalGenerator::Integrate(int l, int m, int n, bool teMode, bool eField)
    {
    	std::vector<double> aField;
    	double r, theta, z = 0.;
    	double dR = fR/fnPixels;
    	double dZ = fL/fnPixels;
    	double dTheta = 2.*LMCConst::Pi()/fnPixels;
    	double tVolume = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<fnPixels; i++)
    		for (unsigned j=0; j<fnPixels; j++)
    			for (unsigned k=0; k<fnPixels; k++)
    			{
    	    		r = (double)i*dR;
    	    		theta = (double)j*dTheta;
    	    		z = (double)k*dZ;
    	    		if (teMode)
    	    		{
    	    			if (eField)
    	    			{
    	    		    	aField = TE_E(l, m, n, r, theta, z);
    	    			}
    	    			else
    	    			{
    	    				aField = TE_H(l, m, n, r, theta, z);
    	    			}
    	    		}
    	    		else
    	    		{
    	    			if (eField)
    	    			{
    	    				aField = TM_E(l, m, n, r, theta, z);
    	    			}
    	    			else
    	    			{
    	    				aField = TM_H(l, m, n, r, theta, z);
    	    			}
    	    		}

    	    		double aFieldMagSq = 0.;
    	    		auto it = aField.begin();
    	    		while (it != aField.end())
    	    		{
		    			if (!isnan(*it))
		    				aFieldMagSq += (*it)*(*it);
    	    			*it++;
    	    		}

    				tIntegral += aFieldMagSq*r*dR*dTheta*dZ;
//    		    	tVolume += r*dR*dTheta*dZ;  // sanity check volume integral.
    			}
//    	printf("tVolume is %g\n", tVolume); getchar();
    	return tIntegral;
    }

    std::vector<double> CavitySignalGenerator::TE_E(int l, int m, int n, double r, double theta, double z) const
    {
    	// from Pozar
    	std::vector<double> TE_E;
    	double x_lm = fInterface->fBesselNKPrimeZeros[l][m];
    	double k1 = x_lm / fR;
    	double k3 = n * LMCConst::Pi() / fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double omega = LMCConst::C()*k;
    	double k0 = omega/LMCConst::C()*sqrt(LMCConst::MuNull()*LMCConst::EpsNull());
    	double eta = LMCConst::MuNull()*omega/LMCConst::C()/k0;  // Jackson 8.32
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double tEr = -l * k/k1 * eta * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
    	double tEtheta = -k/k1 * eta * jPrime * cos(l*theta) * sin(k3*z);
    	TE_E.push_back(tEr);
    	TE_E.push_back(tEtheta);
        return TE_E;
    }

    std::vector<double> CavitySignalGenerator::TE_H(int l, int m, int n, double r, double theta, double z) const
    {
    	// from Pozar
    	std::vector<double> TE_H;
    	double x_lm = fInterface->fBesselNKPrimeZeros[l][m];
    	double k1 = x_lm / fR;
    	double k3 = n * LMCConst::Pi() / fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
    	TE_H.push_back(-k3/k1 * jPrime * cos(l*theta) * cos(k3*z));
    	TE_H.push_back(-l*k3/k1 * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z));
    	TE_H.push_back(boost::math::cyl_bessel_j(l, k1*r) * cos(l*theta) * sin(k3*z));
        return TE_H;
    }

    std::vector<double> CavitySignalGenerator::TM_E(int l, int m, int n, double r, double theta, double z) const
    {
    	// from Pozar
    	std::vector<double> TM_E;
    	double x_lm = fInterface->fBesselNKZeros[l][m];
    	double k1 = x_lm / fR;
    	double k3 = n * LMCConst::Pi() / fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double omega = LMCConst::C()*k;
    	double k0 = omega/LMCConst::C()*sqrt(LMCConst::MuNull()*LMCConst::EpsNull());
    	double eta = LMCConst::C()/LMCConst::EpsNull()/omega*k0;  // Jackson 8.32
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
    	TM_E.push_back(-k3/k1 * eta * jPrime * cos(l*theta) * cos(k3*z));
    	TM_E.push_back(-l*k3/k1 * eta * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z));
    	TM_E.push_back(eta * boost::math::cyl_bessel_j(l, k1*r) * cos(l*theta) * sin(k3*z));
        return TM_E;
    }

    std::vector<double> CavitySignalGenerator::TM_H(int l, int m, int n, double r, double theta, double z) const
    {
    	// from Pozar
    	std::vector<double> TM_H;
    	double x_lm = fInterface->fBesselNKZeros[l][m];
    	double k1 = x_lm / fR;
    	double k3 = n * LMCConst::Pi() / fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double tHr = -l * k/k1  * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
    	double tHtheta = -k/k1 * jPrime * cos(l*theta) * sin(k3*z);
    	TM_H.push_back(tHr);
    	TM_H.push_back(tHtheta);
        return TM_H;
    }



    void CavitySignalGenerator::ReadBesselZeroes(std::string filename, std::vector<std::vector<double> > &data)
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

    void CavitySignalGenerator::CheckNormalization()
    {
    	printf("\\epsilon\\int{|E_xlm|^2 dV} = \\mu\\int{|H_xlm|^2 dV} ?\n\n");
    	for (int l=0; l<3; l++)
    		for (int m=1; m<3; m++)
    			for (int n=1; n<3; n++)
    			{
    		    	printf("TE%d%d%d E %g H %g\n", l, m, n, LMCConst::EpsNull()*Integrate(l,m,n,1,1), LMCConst::MuNull()*Integrate(l,m,n,1,0));
    			}
    	for (int l=0; l<3; l++)
    		for (int m=1; m<3; m++)
    			for (int n=1; n<3; n++)
    			{
    		    	printf("TM%d%d%d E %g H %g\n", l, m, n, LMCConst::EpsNull()*Integrate(l,m,n,0,1), LMCConst::MuNull()*Integrate(l,m,n,0,0));
    			}
    }



    void CavitySignalGenerator::PrintModeMaps()
    {
    	char buffer[60];

    	for (int l=0; l<5; l++)
    		for (int m=1; m<5; m++)
    			for (int n=1; n<5; n++)
    			{
    				printf("l m n is %d %d %d\n", l, m, n);
    				double tNormalizationTE_E = pow(Integrate(l,m,n,1,1),0.5);
    				int a = sprintf(buffer, "output/ModeMapTE%d%d%d_E.txt", l, m, n);
    				const char *fpname = buffer;
    				FILE *fpTE_E = fopen(fpname, "w");
    				for (unsigned i=0; i<fnPixels; i++)
    				{
    					double r = (double)i/fnPixels*fR;
    					for (unsigned j=0; j<fnPixels; j++)
    					{
    						double theta = (double)j/fnPixels*2.*LMCConst::Pi();
    						std::vector<double> tTE_E = TE_E(l,m,n,r,theta,0.05);
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

        fInterface->dtFilter = fInterface->fTFReceiverHandler.GetFilterResolution();
    	FieldBuffer aFieldBuffer;
    	fInterface->eCurrentBuffer = aFieldBuffer.InitializeBuffer(1, 1, fInterface->fTFReceiverHandler.GetFilterSize());


        scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
        ReadBesselZeroes((dataDir / "BesselZeros.txt").string(), fInterface->fBesselNKZeros );
        ReadBesselZeroes((dataDir / "BesselPrimeZeros.txt").string(), fInterface->fBesselNKPrimeZeros );
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


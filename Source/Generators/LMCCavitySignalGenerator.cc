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
		fE_Gun( false ),
		fNModes( 2 ),
        gxml_filename("blank.xml"),
        fphiLO(0.),
		fNPreEventSamples( 150000 ),
		fThreadCheckTime(100),
		fKassNeverStarted( false ),
		fSkippedSamples( false ),
		fBypassTF( false ),
		fNormCheck( false ),
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


    bool CavitySignalGenerator::ModeSelect(int l, int m, int n, bool eGun)
    {
    	if (eGun)
    	{
    		if ((l==0)&&(m==1)&&(n==0))
    			return true;
    		else
    			return false;
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
        	for (unsigned m=0; m<nModes; m++)
        	{
            	for (unsigned n=0; n<nModes; n++)
            	{
            		if (TE)
            		{
            			aModeNormFactor[l][m][n] = 1./fInterface->fField->Integrate(l,m,n,1,1);
            		}
            		else
            		{
            			aModeNormFactor[l][m][n] = 1./fInterface->fField->Integrate(l,m,n,0,1);
            		}
            	}
        	}
    	}

    	return aModeNormFactor;
    }


    void CavitySignalGenerator::CheckNormalization()
    {

    	if (!fE_Gun)
    		printf("\n \\int{|E_xlm|^2 dV} = \\mu / \\epsilon \\int{|H_xlm|^2 dV} ?\n\n");
    	else
    		printf("\n |E_mn|^2 dA = 1.0.  |H_mn| can vary.  Index l is not used.\n\n");

    	for (int l=0; l<fNModes; l++)
    	{
    		for (int m=1; m<fNModes; m++)
    		{
    			for (int n=0; n<fNModes; n++)
    			{
    				double normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n] / LMCConst::EpsNull();
    				if (!isnan(normFactor)&&(isfinite(normFactor)))
    				{
    					printf("TE%d%d%d E %.4g H %.4g\n", l, m, n, LMCConst::EpsNull()*fInterface->fField->Integrate(l,m,n,1,1)*normFactor,
        		    		LMCConst::MuNull()*fInterface->fField->Integrate(l,m,n,1,0)*normFactor);
    				}
    				else
    				{
    					printf("TE%d%d%d is undefined.\n", l, m, n);
    				}

    			}
    		}
    	}


    	for (int l=0; l<fNModes; l++)
    	{
    		for (int m=1; m<fNModes; m++)
    		{
    			for (int n=1; n<fNModes; n++)
    			{
    				double normFactor = fInterface->fField->GetNormFactorsTM()[l][m][n] / LMCConst::EpsNull();
    				if (!isnan(normFactor)&&(isfinite(normFactor)))
    				{
    					printf("TM%d%d%d E %.4g H %.4g\n", l, m, n, LMCConst::EpsNull()*fInterface->fField->Integrate(l,m,n,0,1)*normFactor,
    		    			LMCConst::MuNull()*fInterface->fField->Integrate(l,m,n,0,0)*normFactor);
    				}
    				else
    				{
    					printf("TM%d%d%d is undefined.\n", l, m, n);
    				}
    			}
    		}
    	}

    	printf("\nThe modes normalized as above are available for use in the simulation.\n\n");
    }



    void CavitySignalGenerator::PrintModeMaps()
    {
    	char buffer[60];
    	unsigned modeCounter = 0;

    	for (int l=0; l<3; l++)
    		for (int m=1; m<3; m++)
    			for (int n=0; n<3; n++)
    			{
    				printf("l m n is %d %d %d\n", l, m, n);
    				double tNormalizationTE_E = fInterface->fField->GetNormFactorsTE()[l][m][n];
    				int a = sprintf(buffer, "output/ModeMapTE%d%d%d_E.txt", l, m, n);
    				const char *fpname = buffer;
    				FILE *fpTE_E = fopen(fpname, "w");
    				for (unsigned i=0; i<fInterface->fnPixels; i++)
    				{
    					double r = (double)i/fInterface->fnPixels*fInterface->fR;
    					double x = (double)i/fInterface->fnPixels*fInterface->fX - fInterface->fX/2.;
    					for (unsigned j=0; j<fInterface->fnPixels; j++)
    					{
    						double theta = (double)j/fInterface->fnPixels*2.*LMCConst::Pi();
        					double y = (double)j/fInterface->fnPixels*fInterface->fY - fInterface->fY/2.;
    						std::vector<double> tTE_E;
    						if (!fE_Gun)
    						{
    							tTE_E = fInterface->fField->TE_E(l,m,n,r,theta,0.0,fInterface->fField->GetCentralFrequency());
    							fprintf(fpTE_E, "%10.4g %10.4g %10.4g %10.4g\n", r, theta, tTE_E.front()*tNormalizationTE_E, tTE_E.back()*tNormalizationTE_E);
    						}
    						else
    						{
    							tTE_E = fInterface->fField->TE_E(m,n,x,y,fInterface->fField->GetCentralFrequency());
    							fprintf(fpTE_E, "%10.4g %10.4g %10.4g %10.4g\n", x, y, tTE_E.front()*tNormalizationTE_E, tTE_E.back()*tNormalizationTE_E);
    						}

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

        if( aParam.has( "e-gun" ) )
        {

        	fE_Gun = aParam["e-gun"]().as_bool();

            if (fE_Gun)
            {
            	fInterface->fField = new RectangularWaveguide;
            }
            else
            {
            	fInterface->fField = new CylindricalCavity;
            }

        }

        if( aParam.has( "central-frequency" ) )
        {
            fInterface->fField->SetCentralFrequency(2.*LMCConst::Pi()*aParam["central-frequency"]().as_double());
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

        if( aParam.has( "bypass-tf" ) )
        {
        	fBypassTF = aParam["bypass-tf"]().as_bool();
        }

        if( aParam.has( "xml-filename" ) )
        {
            gxml_filename = aParam["xml-filename"]().as_string();
        }


		fPowerCombiner = new CavityModes;
		fPowerCombiner->SetNCavityModes(fNModes);
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

        scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
        ReadBesselZeroes((dataDir / "BesselZeros.txt").string(), 0 );
        ReadBesselZeroes((dataDir / "BesselPrimeZeros.txt").string(), 1 );

        fInterface->fField->SetNormFactorsTE(CalculateNormFactors(fNModes,1));
        fInterface->fField->SetNormFactorsTM(CalculateNormFactors(fNModes,0));
        CheckNormalization();  // E fields integrate to 1.0 for both TE and TM modes.  H fields near 1.0
//        getchar();
//        PrintModeMaps();
//        getchar();

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
    	// Calculate ratio between mode E-field max amplitude at electron and mode H-field max amplitude
    	// at probe.
    	// Ignore 90-degree phase difference between E-field and H-field as it is constant.  This means
    	// that the data stream will have an artificial phase shift of 90 degrees everywhere, which seems
    	// acceptable if we are only sensitive to phase evolution and phase differences between probes.

    	double fcyc = tKassParticleXP[7];

    	std::vector<double> tTE_E_electron = fInterface->fField->TE_E(0,1,1,tKassParticleXP[0],0.,tKassParticleXP[2], fcyc);
    	double thetaProbe = tKassParticleXP[1] + fPowerCombiner->GetCavityProbeTheta()[channelIndex];
    	double zProbe = fPowerCombiner->GetCavityProbeZ()[channelIndex];
    	std::vector<double> tTE_H_probe = fInterface->fField->TE_H(0,1,1,fInterface->fR,thetaProbe,zProbe, fcyc);
        double modeScalingFactor = tTE_H_probe.back() / tTE_E_electron.back();
        //        	printf("modeScalingFactor is %g / %g = %g\n", tTE_H_probe.back(), tTE_E_electron.back(), modeAmplitudeFactor); getchar();

    	return modeScalingFactor;
    }



    double CavitySignalGenerator::GetWaveguideDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized)
    {
    	double tThetaParticle = tKassParticleXP[1];
    	double tEy = aTE_E_normalized.back();
     	double tEmag = fabs(tEy);
    	// fix-me:  should below be instantaneous velocity or center of motion?
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double tVmag = pow(tVx*tVx + tVy*tVy, 0.5);
//    	printf("tEmag is %g, r is %g, and theta is %g\n", tEmag, tKassParticleXP[0], tKassParticleXP[1]); getchar();
    	return fabs(0. + tEy*tVy)/tEmag/tVmag;  // fabs ( unit J \cdot unit E )
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

    std::vector<double> CavitySignalGenerator::GetWaveguideNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP)
     {
    	// The l index is inert in the waveguide.
     	double tX = tKassParticleXP[0] * sin(tKassParticleXP[1]);
     	double tY = tKassParticleXP[0] * cos(tKassParticleXP[1]);
     	double fcyc = tKassParticleXP[7];
     	std::vector<double> tTE_E_electron = fInterface->fField->TE_E(m,n,tX,tY,fcyc);
 		double normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];  // select mode 0,1,1
// 		printf("normFactor is %g and Ey is %g\n", normFactor, tTE_E_electron.back()); getchar();

 		auto it = tTE_E_electron.begin();
 		while (it != tTE_E_electron.end())
 		{
 			if (!isnan(*it))
 				(*it) *= normFactor;
 			*it++;
 		}
     	return tTE_E_electron;  // return normalized field.
     }


    std::vector<double> CavitySignalGenerator::GetCavityNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP)
     {
     	double tR = tKassParticleXP[0];
     	double tZ = tKassParticleXP[2];
     	double fcyc = tKassParticleXP[7];
     	std::vector<double> tTE_E_electron = fInterface->fField->TE_E(l,m,n,tR,0.,tZ, fcyc);

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
				aLocalElementFIRBuffer[0].push_back(cos(orbitPhase));
			}
			else
			{
				aLocalElementFIRBuffer[0].push_back(0.);
			}
			aLocalElementFIRBuffer[0].pop_front();

			*it++;
		}

		if ( !fBypassTF )
		{
			convolution = fInterface->fTFReceiverHandler.ConvolveWithFIRFilter(aLocalElementFIRBuffer[0]);
		}
		else
		{
			convolution = 1.0;
		}

		return LMCConst::Q()*vMag*convolution;

    }


    bool CavitySignalGenerator::DriveMode(Signal* aSignal, int nFilterBinsRequired, double dtFilter, unsigned index)
    {
        const int signalSize = aSignal->TimeSize();
        unsigned sampleIndex = 0;
        const unsigned nChannels = fNChannels;

        //Receiver Properties
        fphiLO += 2. * LMCConst::Pi() * fLO_Frequency * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

    	std::vector<double> tKassParticleXP = fInterface->fTransmitter->ExtractParticleXP(fInterface->fTOld);
        double dotProductFactor = 0.;
        double impedanceFactor = 0.;
        double unitConversion = 1./sqrt(4.*LMCConst::Pi()*LMCConst::EpsNull()); // Gaussian -> S.I. units

    	for (int l=0; l<fNModes; l++)
    	{
    		for (int m=1; m<fNModes; m++)
    		{
    			for (int n=0; n<fNModes; n++)
    			{
    				if (ModeSelect(l, m, n, fE_Gun))
    				{
    			    	std::vector<double> tTE_E_normalized;
    					impedanceFactor = 2. * LMCConst::Pi() * fInterface->fField->Z_TE(l,m,n) / LMCConst::C(); // Jackson Eq. 8.140
    					if (!fE_Gun)
    					{
    						tTE_E_normalized = GetCavityNormalizedModeField(l,m,n,tKassParticleXP);
    						dotProductFactor = GetCavityDotProductFactor(tKassParticleXP, tTE_E_normalized);  // unit velocity \dot unit theta
    					}
    					else
    					{
    						tTE_E_normalized = GetWaveguideNormalizedModeField(l,m,n,tKassParticleXP);
    						dotProductFactor = GetWaveguideDotProductFactor(tKassParticleXP, tTE_E_normalized);  // unit velocity \dot unit theta
    					}
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

    						// TO-DO:  Implement modeScalingFactor to calculate signal amplitude at probe/channel using
    						// H-field and probe inductance.  Probe inductance is presently set to a fake value of 1.0
    						// in LMCPowerCombiner.cc
    						// double modeScalingFactor = GetModeScalingFactor(tKassParticleXP, channelIndex);
    						double modeScalingFactor = 1.0; // override for now.
    						// end TO-DO

    						// This scaling factor includes a 50 ohm impedance that is expected in signal processing, as well
    						// as other factors as defined above.
    						double totalScalingFactor = sqrt(50.) * unitConversion * impedanceFactor * modeScalingFactor * dotProductFactor * modeAmplitude;
    						fPowerCombiner->AddOneModeToCavityProbe(aSignal, tFirSample, fphiLO, totalScalingFactor, fPowerCombiner->GetCavityProbeInductance(), sampleIndex);
    						if (fNormCheck) fPowerCombiner->AddOneSampleToRollingAvg(l, m, n, tFirSample, totalScalingFactor, sampleIndex);
    					}

    				} // ModeSelect
    			} // n
    		} // m
    	} // l

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


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
        fAvgDotProductFactor({0.}),
        fNPreEventSamples( 150000 ),
        fThreadCheckTime(100),
        fKassNeverStarted( false ),
        fAliasedFrequencies( false ),
        fOverrideAliasing( false ),
        fOverrideStepsize( false ),
        fBypassTF( false ),
        fNormCheck( false ),
        fModeMaps( false ),
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
                        //if ((l==0)&&(m==1)&&(n==1))
    			if ((l==1)&&(m==1)&&(n==1))
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

    	if (!fInterface->fE_Gun)
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
    	char bufferE[60];
    	char bufferH[60];
    	unsigned modeCounter = 0;

    	for (int l=0; l<fNModes; l++)
    		for (int m=1; m<fNModes; m++)
    			for (int n=0; n<fNModes; n++)
    			{
    				printf("l m n is %d %d %d\n", l, m, n);
    				double normFactor = 1.0;
    				int a = 0;
    				if (fTE)
    				{
    					normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];
        				a = sprintf(bufferE, "output/ModeMapTE%d%d%d_E.txt", l, m, n);
        				a = sprintf(bufferH, "output/ModeMapTE%d%d%d_H.txt", l, m, n);
    				}
    				else
    				{
    					normFactor = fInterface->fField->GetNormFactorsTM()[l][m][n];
        				a = sprintf(bufferE, "output/ModeMapTM%d%d%d_E.txt", l, m, n);
        				a = sprintf(bufferH, "output/ModeMapTM%d%d%d_H.txt", l, m, n);
    				}
    				const char *fpnameE = bufferE;
    				FILE *fp_E = fopen(fpnameE, "w");
    				const char *fpnameH = bufferH;
    				FILE *fp_H = fopen(fpnameH, "w");
    				for (unsigned i=0; i<fInterface->fField->GetNPixels()+1; i++)
    				{
    					double r = (double)i/fInterface->fField->GetNPixels()*fInterface->fField->GetDimR();
    					double x = (double)i/fInterface->fField->GetNPixels()*fInterface->fField->GetDimX() - fInterface->fField->GetDimX()/2.;
    					for (unsigned j=0; j<fInterface->fField->GetNPixels()+1; j++)
    					{
    						double theta = (double)j/fInterface->fField->GetNPixels()*2.*LMCConst::Pi();
        					double y = (double)j/fInterface->fField->GetNPixels()*fInterface->fField->GetDimY() - fInterface->fField->GetDimY()/2.;
    						std::vector<double> tE;
    						std::vector<double> tH;
    						if (!fInterface->fE_Gun)
    						{
    							if (fTE)
    							{
    								tE = fInterface->fField->TE_E(l,m,n,r,theta,0.0,0);
    								tH = fInterface->fField->TE_H(l,m,n,r,theta,0.0,0);
    							}
    							else
    							{
    								tE = fInterface->fField->TM_E(l,m,n,r,theta,0.0,0);
    								tH = fInterface->fField->TM_H(l,m,n,r,theta,0.0,0);
    							}
    							fprintf(fp_E, "%10.4g %10.4g %10.4g %10.4g\n", r, theta, tE.front()*normFactor, tE.back()*normFactor);
    							fprintf(fp_H, "%10.4g %10.4g %10.4g %10.4g\n", r, theta, tH.front()*normFactor, tH.back()*normFactor);
    						}
    						else
    						{
    							tE = fInterface->fField->TE_E(m,n,x,y,fInterface->fField->GetCentralFrequency());
    							fprintf(fp_E, "%10.4g %10.4g %10.4g %10.4g\n", x, y, tE.front()*normFactor, tE.back()*normFactor);
    						}
    					}
    				}
    				fclose (fp_E);
    				fclose (fp_H);
    				modeCounter += 1;
    			}
    	printf("\nMode map files have been generated; press RETURN to continue, or Cntrl-C to quit.\n");
    	getchar();
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

        if (!fInterface->fField->Configure(aParam))
        {
        	LERROR(lmclog,"Error configuring LMCField.");
        	exit(-1);
        	return false;
        }

        if( aParam.has( "cavity-radius" ) )
        {
            fInterface->fField->SetDimR( aParam["cavity-radius"]().as_double() );
        }

        if( aParam.has( "cavity-length" ) )
        {
            fInterface->fField->SetDimL( aParam["cavity-length"]().as_double() );
        }

        if( aParam.has( "center-to-short" ) ) // for use in e-gun
        {
            fInterface->fCENTER_TO_SHORT = aParam["center-to-short"]().as_double();
        }

        if( aParam.has( "center-to-antenna" ) ) // for use in e-gun
        {
            fInterface->fCENTER_TO_ANTENNA = aParam["center-to-antenna"]().as_double();
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
        if( aParam.has( "mode-maps" ) )
        {
        	fModeMaps = aParam["mode-maps"]().as_bool();
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


        fdtFilter = fTFReceiverHandler->GetFilterResolution();

        scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
        ReadBesselZeroes((dataDir / "BesselZeros.txt").string(), 0 );
        ReadBesselZeroes((dataDir / "BesselPrimeZeros.txt").string(), 1 );

        fInterface->fField->SetNormFactorsTE(CalculateNormFactors(fNModes,1));
        fInterface->fField->SetNormFactorsTM(CalculateNormFactors(fNModes,0));
        CheckNormalization();  // E fields integrate to 1.0 for both TE and TM modes.
        if (fModeMaps) PrintModeMaps();  // Write to text files for plotting mode maps.

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
	std::vector<double> excitationAmplitude = {0.};
	std::vector<double> tEFieldAtProbe = {0.};
        std::vector<double> dopplerFrequency;


    	for (int l=0; l<fNModes; l++)
    	{
		fAvgDotProductFactor.resize(l+1,0.);
		tEFieldAtProbe.resize(l+1,0.);
		excitationAmplitude.resize(l+1,0.);
    		for (int m=1; m<fNModes; m++)
    		{
    			for (int n=0; n<fNModes; n++)
    			{
    				if (ModeSelect(l, m, n, fInterface->fE_Gun))
    				{
					std::vector<std::vector<double>> tE_normalized;
    					tE_normalized = fInterface->fField->GetNormalizedModeFields(l,m,n,tKassParticleXP);
    					double cavityFIRSample = fFieldCalculator->GetCavityFIRSample(tKassParticleXP, fBypassTF).first;
    					dopplerFrequency = fInterface->fField->GetDopplerFrequency(l, m, n, tKassParticleXP);
					std::vector<double> modeAmplitude(l+1,0.0);
					for(int j=0; j<=l; j++)
					{
						fAvgDotProductFactor[j] = fInterface->fField->GetDotProductFactor(tKassParticleXP, tE_normalized[j], fIntermediateFile);
	                                        //fAvgDotProductFactor[j] = 1. / ( tThisEventNSamples + 1 ) * ( fAvgDotProductFactor[j] * tThisEventNSamples + fInterface->fField->GetDotProductFactor(tKassParticleXP, tE_normalized[j], fIntermediateFile) );  // unit velocity \dot unit theta
						double modeSign = tE_normalized[j].front() / fabs(tE_normalized[j].front());
						//modeAmplitude[j] = modeSign*sqrt( tE_normalized[j].back()*tE_normalized[j].back() + tE_normalized[j].front()*tE_normalized[j].front());
						modeAmplitude[j] = tE_normalized[j].back();

    						if (!fInterface->fE_Gun)
    						{
    							// sqrt(4PIeps0) for Kass current si->cgs, sqrt(4PIeps0) for Jackson A_lambda coefficient cgs->si
    							unitConversion = 1. / LMCConst::FourPiEps(); // see comment ^
    							// Calculate propagating E-field with J \dot E.  cavityFIRSample units are [current]*[unitless].
							//std::cout << "Phi, DotProductFactor, ModeAmplitude: " << tKassParticleXP[1] << " " << fAvgDotProductFactor[j] <<	" " << modeAmplitude[j] << std::endl;
							excitationAmplitude[j] = fAvgDotProductFactor[j] * modeAmplitude[j] * cavityFIRSample * fInterface->fField->Z_TE(l,m,n,tKassParticleXP[7]) * 2. * LMCConst::Pi() / LMCConst::C() / 1.e2;
    							std::vector<double> tProbeLocation = {fInterface->fField->GetDimR()*fPowerCombiner->GetCavityProbeRFrac(), fPowerCombiner->GetCavityProbePhi(), fPowerCombiner->GetCavityProbeZ()};
    							tEFieldAtProbe[j] = fInterface->fField->GetNormalizedModeFields(l,m,n,tProbeLocation)[j].back();
    						}
    						else
    						{
    							// sqrt(4PIeps0) for Kass current si->cgs, sqrt(4PIeps0) for Jackson A_lambda coefficient cgs->si
    							unitConversion = 1. / LMCConst::FourPiEps(); // see comment ^
    							// Calculate propagating E-field with J \dot E and integrated Poynting vector.  cavityFIRSample units are [current]*[ohms].
    							excitationAmplitude[j] = fAvgDotProductFactor[j] * modeAmplitude[j] * ScaleEPoyntingVector(tKassParticleXP[7]) * cavityFIRSample * 2. * LMCConst::Pi() / LMCConst::C() / 1.e2;

    							// Optional cross-check:  Use direct Kassiopeia power budget.  Assume x_electron = 0.
    							if (fUseDirectKassPower)
    							{
    								// override one-way signal amplitude with direct Kass power:
    								unitConversion = 1.0;  // Kass power is already in Watts.
        							excitationAmplitude[j] = fAvgDotProductFactor[j]*sqrt(tKassParticleXP[8]/2.);  // sqrt( modeFraction*LarmorPower/2 )
    							}

    						}
					}
    					for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)  // one channel per probe.
    					{
    						sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample

    						// This scaling factor includes a 50 ohm impedance that applied in signal processing, as well
    						// as other factors as defined above, e.g. 1/4PiEps0 if converting to/from c.g.s amplitudes.
    						double totalScalingFactor = sqrt(50.) * unitConversion;
						//std::cout << std::endl<< "excitationAmplitude: " << excitationAmplitude[0] << ", " << excitationAmplitude[1] << std::endl;
						//std::cout << "tEFieldAtProbe: " << tEFieldAtProbe[0] << ", " << tEFieldAtProbe[1] << std::endl;
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


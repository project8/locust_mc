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
		fModeMaps( false ),
		fTE( true ),
		fIntermediateFile( false ),
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
        				a = sprintf(buffer, "output/ModeMapTE%d%d%d_E.txt", l, m, n);
    				}
    				else
    				{
    					normFactor = fInterface->fField->GetNormFactorsTM()[l][m][n];
        				int a = sprintf(buffer, "output/ModeMapTM%d%d%d_E.txt", l, m, n);
    				}
    				const char *fpname = buffer;
    				FILE *fp_E = fopen(fpname, "w");
    				for (unsigned i=0; i<fInterface->fField->GetNPixels()+1; i++)
    				{
    					double r = (double)i/fInterface->fField->GetNPixels()*fInterface->fR;
    					double x = (double)i/fInterface->fField->GetNPixels()*fInterface->fX - fInterface->fX/2.;
    					for (unsigned j=0; j<fInterface->fField->GetNPixels()+1; j++)
    					{
    						double theta = (double)j/fInterface->fField->GetNPixels()*2.*LMCConst::Pi();
        					double y = (double)j/fInterface->fField->GetNPixels()*fInterface->fY - fInterface->fY/2.;
    						std::vector<double> tE;
    						if (!fE_Gun)
    						{
    							if (fTE)
    							{
    								tE = fInterface->fField->TE_E(l,m,n,r,theta,0.0,fInterface->fField->GetCentralFrequency());
    							}
    							else
    							{
    								tE = fInterface->fField->TM_E(l,m,n,r,theta,0.0,fInterface->fField->GetCentralFrequency());
    							}
    							fprintf(fp_E, "%10.4g %10.4g %10.4g %10.4g\n", r, theta, tE.front()*normFactor, tE.back()*normFactor);
    						}
    						else
    						{
    							tE = fInterface->fField->TE_E(m,n,x,y,fInterface->fField->GetCentralFrequency());
    							fprintf(fp_E, "%10.4g %10.4g %10.4g %10.4g\n", x, y, tE.front()*normFactor, tE.back()*normFactor);
    						}

    					}
    				}
    				fclose (fp_E);
    				modeCounter += 1;
    			}
    	printf("\nMode map files have been generated; press RETURN to continue, or Cntrl-C to quit.\n");
    	getchar();
    }



    bool CavitySignalGenerator::Configure( const scarab::param_node& aParam )
    {

    	if(!fInterface->fTFReceiverHandler.Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver FIRHandler class");
    	}


// To-do:  Move this into EquivalentCircuit::Configure()
//-------------New implementation of fEquivalentCircuit for cavity parameterization-----
	if ( (aParam.has( "equivalentR" )) or (aParam.has( "equivalentL" )) or (aParam.has( "equivalentL" )) ) { //Only initiate the configurable TF if at least one parameter is specified in the config file
		//Default RLC parameters if not defined in config file, values (mostly arbitratily) based on external script from P. Slocum for 1 GHz cavity
		printf("Entering RLC Config Loop\n");
		double equivalentR = 1.;
                double equivalentL = 0.159e-6;
                double equivalentC = 0.159e-12;
		int TFBins = 4000;
		double FreqRangeCenter = 1.0e9;

		//Update any parameters defined in the config file
                if( aParam.has( "equivalentR" ) )
                {
			//printf("Retrieved R value: %e\n", aParam["equivalentR"]().as_double());
                        equivalentR = aParam["equivalentR"]().as_double();
                }
                if( aParam.has( "equivalentL" ) )
                {
                        //printf("Retrieved L value: %e\n", aParam["equivalentL"]().as_double());
                        equivalentL = aParam["equivalentL"]().as_double();
                }
                if( aParam.has( "equivalentC" ) )
                {
                        //printf("Retrieved C value: %e\n", aParam["equivalentC"]().as_double());
                        equivalentC = aParam["equivalentC"]().as_double();
                }
		if( aParam.has( "TFBins" ) )
		{
			//printf("Retrieved TFBins: %d\n", aParam["TFBins"]().as_int());
			TFBins = aParam["TFBins"]().as_int();
		}
		if( aParam.has( "FreqRangeCenter" ) )
		{
			//printf("Retrieved FreqRangeCenter: %f\n", aParam["FreqRangeCenter"]().as_double());
                        FreqRangeCenter = aParam["FreqRangeCenter"]().as_double();
		}

		//Call EquivalentCircuit Class and use to analytically generate a Transfer function
                fEquivalentCircuit = new EquivalentCircuit();
                fEquivalentCircuit->GenerateTransferFunction(equivalentR,equivalentL,equivalentC,TFBins,FreqRangeCenter); //R,L,C inputs

		if(!fInterface->fTFReceiverHandler.ConvertAnalyticTFtoFIR(fEquivalentCircuit->initialFreq,fEquivalentCircuit->tfArray))
		{
		return false;
		}
	}
	else{
        	if(!fInterface->fTFReceiverHandler.ReadHFSSFile())
        	{
		printf("Using externally built Transfer Function via file \n");
            	return false;
        	}
	}
//--------------------------------------------------------------------------------------

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
        CheckNormalization();  // E fields integrate to 1.0 for both TE and TM modes.
        if (fModeMaps) PrintModeMaps();  // Write to text files for plotting mode maps.

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



    double CavitySignalGenerator::GetWaveguideDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized)
    {
    	double tThetaParticle = tKassParticleXP[1];
    	double tEy = aTE_E_normalized.back();
     	double tEmag = fabs(tEy);
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double tVmag = pow(tVx*tVx + tVy*tVy, 0.5);
    	double unitJdotE = fabs(0. + tEy*tVy)/tEmag/tVmag;


    	//  Write trajectory points, dot product, and E-field mag to file for debugging etc.
    	if (fIntermediateFile)
    	{
        	char buffer[60];
    		int a = sprintf(buffer, "output/dotProducts.txt");
    		const char *fpname = buffer;
    		FILE *fp = fopen(fpname, "a");
    		fprintf(fp, "%g %g %g %g\n", tKassParticleXP[0], tKassParticleXP[1], unitJdotE, tEmag);
    		fclose(fp);

    		printf("unitJdotE is %g, r*cos(theta) is %g, r is %g, and theta is %g, eMag is %g\n",
    			unitJdotE, tKassParticleXP[0]*cos(tKassParticleXP[1]), tKassParticleXP[0], tKassParticleXP[1], tEmag); getchar();
    	}

    	return unitJdotE;
    }

    double CavitySignalGenerator::GetCavityDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> anE_normalized)
    {
    	double tThetaParticle = tKassParticleXP[1];
    	double tEtheta = 0.;
    	double tEr = 0.;
    	if (!isnan(anE_normalized.back()))
    	{
    		tEtheta = anE_normalized.back();
    	}
    	if (!isnan(anE_normalized.front()))
    	{
    		tEr = anE_normalized.front();
    	}
    	double tEx = -sin(tThetaParticle) * tEtheta + cos(tThetaParticle) * tEr;
    	double tEy = cos(tThetaParticle) * tEtheta + sin(tThetaParticle) * tEr;
    	double tEmag = pow(tEtheta*tEtheta + tEr*tEr, 0.5);
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double tVmag = pow(tVx*tVx + tVy*tVy, 0.5);
    	double unitJdotE = fabs(tEx*tVx + tEy*tVy)/tEmag/tVmag;


    	//  Write trajectory points, dot product, and E-field mag to file for debugging etc.
    	if (fIntermediateFile)
    	{
        	char buffer[60];
    		int a = sprintf(buffer, "output/dotProducts.txt");
    		const char *fpname = buffer;
    		FILE *fp = fopen(fpname, "a");
    		fprintf(fp, "%g %g %g %g\n", tKassParticleXP[0], tKassParticleXP[1], unitJdotE, tEmag);
    		fclose(fp);

    		printf("unitJdotE is %g, r*cos(theta) is %g, r is %g, and theta is %g, eMag is %g\n",
    			unitJdotE, tKassParticleXP[0]*cos(tKassParticleXP[1]), tKassParticleXP[0], tKassParticleXP[1], tEmag); getchar();
    	}

    	return unitJdotE;
    }

    std::vector<double> CavitySignalGenerator::GetWaveguideNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP)
     {
    	// The l index is inert in the waveguide.
     	double tX = tKassParticleXP[0] * cos(tKassParticleXP[1]);
     	double tY = tKassParticleXP[0] * sin(tKassParticleXP[1]);
     	double fcyc = tKassParticleXP[7];
     	std::vector<double> tTE_E_electron = fInterface->fField->TE_E(m,n,tX,tY,fcyc);
 		double normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];

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
     	std::vector<double> tE_electron;
     	double normFactor = 0.;

     	if (fTE)
     	{
     		tE_electron = fInterface->fField->TE_E(l,m,n,tR,0.,tZ, fcyc);
     		normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];
     	}
     	else
     	{
     		tE_electron = fInterface->fField->TM_E(l,m,n,tR,0.,tZ, fcyc);
     		normFactor = fInterface->fField->GetNormFactorsTM()[l][m][n];
     	}


 		auto it = tE_electron.begin();

 		while (it != tE_electron.end())
 		{
 			if (!isnan(*it))
 				(*it) *= normFactor;
 			*it++;
 		}

     	return tE_electron;  // return normalized field.
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

		if (( !fBypassTF )&&(!fE_Gun))
		{
			convolution = fInterface->fTFReceiverHandler.ConvolveWithFIRFilter(aLocalElementFIRBuffer[0]);
		}
		else
		{
			convolution = 1.0;
		}

		return LMCConst::Q()*vMag*convolution;

    }

    double CavitySignalGenerator::ScaleEPoyntingVector(double fcyc)
    {
    	// This function calculates the coefficients of the Poynting vector integral
    	// in the TE10 mode in WR42.  It then returns the sqrt of the half of the propagating
    	// power that is moving toward the antenna.
    	// After Pozar p. 114:
    	double k = fcyc / LMCConst::C();
    	double k1 = LMCConst::Pi() / fInterface->fX;
    	double beta = sqrt( k*k - k1*k1 );
    	double areaIntegral = fcyc * LMCConst::MuNull() * pow(fInterface->fX,3.) * fInterface->fY * beta / 4. / LMCConst::Pi() / LMCConst::Pi();
    	// sqrt of propagating power gives amplitude of E
    	return sqrt(areaIntegral/2.);  // areaIntegral/2. propagates power/2.
    }

    bool CavitySignalGenerator::DriveMode(Signal* aSignal, int nFilterBinsRequired, double dtFilter, unsigned index)
    {
        const int signalSize = aSignal->TimeSize();
        unsigned sampleIndex = 0;
        const unsigned nChannels = fNChannels;

        //Receiver Properties
        double dt = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        fphiLO += 2. * LMCConst::Pi() * fLO_Frequency * dt;

    	std::vector<double> tKassParticleXP = fInterface->fTransmitter->ExtractParticleXP(fInterface->fTOld);
        double dotProductFactor = 0.;
        double unitConversion = 1.;
        double excitationAmplitude = 0.;

    	for (int l=0; l<fNModes; l++)
    	{
    		for (int m=1; m<fNModes; m++)
    		{
    			for (int n=0; n<fNModes; n++)
    			{
    				if (ModeSelect(l, m, n, fE_Gun))
    				{
    			    	std::vector<double> tE_normalized;
    					if (!fE_Gun)
    					{
    						unitConversion = 1.;  // mks units in Collin amplitudes.
    						tE_normalized = GetCavityNormalizedModeField(l,m,n,tKassParticleXP);
    						dotProductFactor = GetCavityDotProductFactor(tKassParticleXP, tE_normalized);  // unit velocity \dot unit theta
    					}
    					else
    					{
    				        // sqrt(4PIeps0) for Kass current si->cgs, sqrt(4PIeps0) for A_lambda coefficient cgs->si
    				        unitConversion = 1. / LMCConst::FourPiEps(); // see comment ^
    						tE_normalized = GetWaveguideNormalizedModeField(l,m,n,tKassParticleXP);
    						dotProductFactor = GetWaveguideDotProductFactor(tKassParticleXP, tE_normalized);  // unit velocity \dot unit theta
    					}

    					// fix-me:  We may need more precise nFilterBinsRequired.
    					std::vector<std::deque<double>> tLocalFIRfrequencyBuffer = fInterface->FIRfrequencyBufferCopy;  // copy from Kass buffer.
    					std::vector<std::deque<double>> tLocalElementFIRBuffer = fInterface->ElementFIRBufferCopy;

    					double modeAmplitude = 0.;
    					if ( (!isnan(tE_normalized.back())) && (!isnan(tE_normalized.front())) )
    					{
    						modeAmplitude = pow(tE_normalized.back()*tE_normalized.back() + tE_normalized.front()*tE_normalized.front(), 0.5);  // normalized E at electron
    					}
    					else if ( !isnan(tE_normalized.back()) )
    					{
    						modeAmplitude = tE_normalized.back();
    					}
    			    	double tDopplerFrequency = fInterface->fField->GetDopplerFrequency(l, m, n, tKassParticleXP);
    					double cavityFIRSample = GetCavityFIRSample(tKassParticleXP, tLocalFIRfrequencyBuffer, tLocalElementFIRBuffer, fInterface->nFilterBinsRequired, fInterface->dtFilter);

    					if (!fE_Gun)
    					{
    						double collinAmplitude = 0.;
    						if (fTE)
    						{
    							collinAmplitude = fInterface->fField->Z_TE(l,m,n,tKassParticleXP[7]);
    						}
    						else
    						{
    							collinAmplitude = fInterface->fField->Z_TM(l,m,n,tKassParticleXP[7]);
    						}
    						excitationAmplitude = modeAmplitude * dotProductFactor * collinAmplitude * cavityFIRSample;
    					}
    					else
    					{
    						// Calculate propagating E-field with J \dot E and integrated Poynting vector:
    						excitationAmplitude = modeAmplitude * dotProductFactor * ScaleEPoyntingVector(tKassParticleXP[7]) *
    							fInterface->fField->Z_TE(l,m,n,tKassParticleXP[7]) * cavityFIRSample;
    						// tExcitationAmplitude = sqrt(tKassParticleXP[8]/2.);  // optional:  unitConversion =1., sqrt( Larmor power / 2 )
    					}
    					std::vector<std::deque<double>>().swap(tLocalFIRfrequencyBuffer);  // release memory
    					std::vector<std::deque<double>>().swap(tLocalElementFIRBuffer);

    					for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)  // one channel per probe.
    					{
    						sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample

    						// This scaling factor includes a 50 ohm impedance that applied in signal processing, as well
    						// as other factors as defined above, e.g. 1/4PiEps0 if converting to/from c.g.s amplitudes.
    						double totalScalingFactor = sqrt(50.) * unitConversion;
    						fPowerCombiner->AddOneModeToCavityProbe(aSignal, excitationAmplitude, tDopplerFrequency, dt, fphiLO, totalScalingFactor, fPowerCombiner->GetCavityProbeInductance(), sampleIndex);
    						if (fNormCheck) fPowerCombiner->AddOneSampleToRollingAvg(l, m, n, excitationAmplitude, sampleIndex);
    					}

    				} // ModeSelect
    			} // n
    		} // m
    	} // l

    	fInterface->fTOld += dt;

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


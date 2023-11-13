/*
 * LMCCylindricalCavity.cc
 *
 *  Created on: Jun 9, 2021
 *      Author: pslocum
 */

#include "LMCCylindricalCavity.hh"


namespace locust
{

    LOGGER( lmclog, "CylindricalCavity" );
    CylindricalCavity::CylindricalCavity():
    	fProbeGain( {1., 1.}),
		fCavityProbeZ( {0., 0.} ),
		fCavityProbeRFrac( {0.5, 0.5} ),
		fCavityProbeTheta( {0.0, 0.0} )
		{}

    CylindricalCavity::~CylindricalCavity() {}



    bool CylindricalCavity::Configure( const scarab::param_node& aParam)
    {

    	if( !Field::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring Field class from CylindricalCavity subclass");
    		return false;
    	}

        if( aParam.has( "cavity-radius" ) )
        {
            SetDimR( aParam["cavity-radius"]().as_double() );
        }

        if( aParam.has( "cavity-length" ) )
        {
        	SetDimL( aParam["cavity-length"]().as_double() );
        }

        if ( aParam.has( "cavity-probe-gain0" ) )
    	{
    		SetCavityProbeGain(aParam["cavity-probe-gain0"]().as_double(), 0);
    	}

        if ( aParam.has( "cavity-probe-gain1" ) )
    	{
    		SetCavityProbeGain(aParam["cavity-probe-gain1"]().as_double(), 1);
    	}

        if ( aParam.has( "cavity-probe-z0" ) )
    	{
    		SetCavityProbeZ(aParam["cavity-probe-z0"]().as_double(), 0);
    	}

        if ( aParam.has( "cavity-probe-z1" ) )
    	{
    		SetCavityProbeZ(aParam["cavity-probe-z1"]().as_double(), 1);
    	}

        if ( aParam.has( "cavity-probe-r-fraction0" ) )
    	{
    		SetCavityProbeRFrac(aParam["cavity-probe-r-fraction0"]().as_double(), 0);
    	}

        if ( aParam.has( "cavity-probe-r-fraction1" ) )
    	{
    		SetCavityProbeRFrac(aParam["cavity-probe-r-fraction1"]().as_double(), 1);
    	}

     	if ( aParam.has( "cavity-probe-theta0" ) )
     	{
     		SetCavityProbeTheta(aParam["cavity-probe-theta0"]().as_double(), 0);
     	}

     	if ( aParam.has( "cavity-probe-theta1" ) )
     	{
     		SetCavityProbeTheta(aParam["cavity-probe-theta1"]().as_double(), 1);
     	}

     	if( aParam.has( "modemap-filename" ) )
     	{
     		fFieldCore = new ModeMapCylindricalCavity();
     		if (!fFieldCore->ReadModeMapTE_E(aParam["modemap-filename"]().as_string()))
     		{
     			LERROR(lmclog,"There was a problem uploading the mode map.");
     			exit(-1);
     		}
     	}
     	else
     	{
     		fFieldCore = new PozarCylindricalCavity();
     		scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
     		fFieldCore->ReadBesselZeroes((dataDir / "BesselZeros.txt").string(), 0 );
     		fFieldCore->ReadBesselZeroes((dataDir / "BesselPrimeZeros.txt").string(), 1 );
     	}

        SetNormFactorsTE(CalculateNormFactors(GetNModes(),1));
        SetNormFactorsTM(CalculateNormFactors(GetNModes(),0));

        CheckNormalization(GetNModes());  // E fields integrate to 1.0 for both TE and TM modes.

        if( PlotModeMaps() )
        {
        	double zSlice = 0.0;
        	double thetaSlice = 0.0;
        	if (aParam.has( "map-z-slice" )) zSlice = aParam["map-z-slice"]().as_double();
        	if (aParam.has( "map-theta-slice" )) thetaSlice = aParam["map-theta-slice"]().as_double();
        	LPROG( lmclog, "If ROOT is available, plotting mode maps to file output/ModemapOutput*.root... " );
        	PrintModeMaps(GetNModes(), zSlice, thetaSlice);
        }

    	return true;
    }

    std::vector<std::vector<std::vector<double>>> CylindricalCavity::CalculateNormFactors(int nModes, bool bTE)
    {

        LPROG(lmclog, "Calculating mode normalization factors ... " );

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
            		if (bTE)
            		{
            			aModeNormFactor[l][m][n] = 1./pow(Integrate(l,m,n,1,1),0.5);
            		}
            		else
            		{
            			aModeNormFactor[l][m][n] = 1./pow(Integrate(l,m,n,0,1),0.5);
            		}
            	}
        	}
    	}

    	return aModeNormFactor;
    }



    double CylindricalCavity::Integrate(int l, int m, int n, bool teMode, bool eField)
    {

    	std::vector<double> aField;
    	double r, theta, zPozar, zKass = 0.;
    	double dR = GetDimR()/GetNPixels();
    	double dZ = GetDimL()/GetNPixels();
    	double dTheta = 2.*LMCConst::Pi()/GetNPixels();
    	double tVolume = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<GetNPixels(); i++)
    	{
    		for (unsigned j=0; j<GetNPixels(); j++)
    		{
    			for (unsigned k=0; k<GetNPixels(); k++)
    			{

    	    		r = (double)i*dR;
    	    		theta = (double)j*dTheta;
    	    		zPozar = (double)k*dZ;
    	    		zKass = zPozar - GetDimL()/2.;

    	    		if (teMode)
    	    		{
    	    			if (eField)
    	    			{
    	    		    	aField = fFieldCore->TE_E(GetDimR(), GetDimL(), l, m, n, r, theta, zKass,0);
    	    			}
    	    			else
    	    			{
    	    				aField = fFieldCore->TE_H(GetDimR(), GetDimL(), l, m, n, r, theta, zKass,0);
    	    			}
    	    		}
    	    		else
    	    		{
    	    			if (eField)
    	    			{
    	    				aField = fFieldCore->TM_E(GetDimR(), GetDimL(), l, m, n, r, theta, zKass,0);
    	    			}
    	    			else
    	    			{
    	    				aField = fFieldCore->TM_H(GetDimR(), GetDimL(), l, m, n, r, theta, zKass,0);
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
    		}
    	}
    	return tIntegral;
    }

    std::vector<double> CylindricalCavity::GetDopplerFrequency(int bTE, int l, int m, int n, std::vector<double> tKassParticleXP)
    {
    	std::vector<double> freqPrime;
    	double vz = tKassParticleXP[5];
	double term1 = 0;
    	if(bTE==1) term1 = fFieldCore->GetBesselNKPrimeZeros(l,m) / GetDimR();
	else term1 = fFieldCore->GetBesselNKZeros(l,m) / GetDimR();
    	double term2 = n * LMCConst::Pi() / GetDimL();
    	double lambda = 1. / pow( 1. / 4. / LMCConst::Pi() / LMCConst::Pi() * ( term1*term1 + term2*term2 ), 0.5);
    	double lambda_c = 0.;
	if(bTE==1)lambda_c = 2 * LMCConst::Pi() * GetDimR() / fFieldCore->GetBesselNKPrimeZeros(l,m);
	else lambda_c = 2 * LMCConst::Pi() * GetDimR() / fFieldCore->GetBesselNKZeros(l,m);
    	double vp = LMCConst::C() / pow( 1. - lambda*lambda/lambda_c/lambda_c, 0.5 );
    	double dopplerShift = 0.;
    	if (vp > 0.) dopplerShift = vz / vp;
		freqPrime.push_back( ( 1. + dopplerShift ) * tKassParticleXP[7] );
    	return freqPrime;
    }


    double CylindricalCavity::Z_TE(int l, int m, int n, double fcyc) const
    {
    	/*
    	double Z_TE = 1.0;
    	double x_lm = fInterface->fBesselNKPrimeZeros[l][m];
    	double k1 = x_lm / GetDimR();
    	double k3 = n * LMCConst::Pi() / GetDimL();
    	double k = pow(k1*k1+k3*k3,0.5);
    	double k0 = fcyc / LMCConst::C();

    	if ( k*k-k0*k0 != 0. )
    	{
    		// Red Jackson wave impedance Eq. 8.32
    		Z_TE = LMCConst::MuNull() * fcyc / LMCConst::C() / k;
    	}
    	*/

    	// Bypass all realistic impedance calculations for now, and return a constant.
    	// This allows us to implement a power model elsewhere in the code without having
    	// a duplicate resonant impedance from this function.
    	double Z_TE = 500.; // ohms
    	return Z_TE;

    }

    double CylindricalCavity::Z_TM(int l, int m, int n, double fcyc) const
    {
    	/*
    	double Z_TM = 1.0;
    	double x_lm = fInterface->fBesselNKZeros[l][m];
    	double k1 = x_lm / GetDimR();
    	double k3 = n * LMCConst::Pi() / GetDimL();
    	double k = pow(k1*k1+k3*k3,0.5);
    	double k0 = fcyc / LMCConst::C();

    	if ( k*k-k0*k0 != 0. )
    	{
    		// Red Jackson waveguide impedance Eq. 8.32
    		Z_TM = 1. / ( LMCConst::EpsNull() * fcyc / LMCConst::C() / k ); // cgs units
    	}
    	*/

    	// See comments in Z_TE, above ^
    	double Z_TM = 500.; //ohms
    	return Z_TM;
    }


    std::vector<double> CylindricalCavity::GetTE_E(int l, int m, int n, double r, double theta, double z, bool includeOtherPols)
    {
    	return fFieldCore->TE_E(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
    }

    std::vector<double> CylindricalCavity::GetTM_E(int l, int m, int n, double r, double theta, double z, bool includeOtherPols)
    {
    	return fFieldCore->TM_E(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
    }

    std::vector<double> CylindricalCavity::GetFieldAtProbe(int l, int m, int n, bool includeOtherPols, std::vector<double> tKassParticleXP, bool teMode)
    {

    	std::vector<double> rProbe;
        rProbe.push_back(GetCavityProbeRFrac()[0] * GetDimR());
        rProbe.push_back(GetCavityProbeRFrac()[1] * GetDimR());

    	std::vector<double> thetaProbe = GetCavityProbeTheta();
    	std::vector<double> zProbe = GetCavityProbeZ();
    	std::vector<double> thetaEffective;

    	if (l>0)
    	{
                //If mode has theta dependence, mode polarization is set by electron location. Probe coupling must be set relative to that angle
                double thetaElectron = tKassParticleXP[1];
                thetaEffective.push_back(thetaProbe[0] - thetaElectron);
                thetaEffective.push_back(thetaProbe[1] - thetaElectron);
    	}
    	else
    	{
    		thetaEffective = thetaProbe;
    	}

    	std::vector<std::vector<double>> tProbeLocation;
    	tProbeLocation.push_back({rProbe[0], thetaEffective[0], zProbe[0]});
    	tProbeLocation.push_back({rProbe[1], thetaEffective[1], zProbe[1]});

    	//Assumes probe couples to E of mode. If mode is polarized, transforms angle to reference frame of electron
    	std::vector<double> tEFieldAtProbe;
    	tEFieldAtProbe.push_back( NormalizedEFieldMag(GetNormalizedModeField(l,m,n,tProbeLocation[0],0,teMode)) );
    	tEFieldAtProbe.push_back( NormalizedEFieldMag(GetNormalizedModeField(l,m,n,tProbeLocation[1],0,teMode)) );

    	return {fProbeGain[0] * tEFieldAtProbe[0], fProbeGain[1] * tEFieldAtProbe[1]};

    }

    std::vector<double> CylindricalCavity::GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP, bool includeOtherPols, bool teMode)
    {
    	double tR = tKassParticleXP[0];
    	double tTheta = tKassParticleXP[1];
    	double tZ = tKassParticleXP[2];
       	std::vector<double> tField;
       	double normFactor;
       	if(teMode)
       	{
       		tField = fFieldCore->TE_E(GetDimR(),GetDimL(),l,m,n,tR,tTheta,tZ,includeOtherPols);
       		normFactor = GetNormFactorsTE()[l][m][n];
       	}
       	else
       	{
       		tField = fFieldCore->TM_E(GetDimR(),GetDimL(),l,m,n,tR,tTheta,tZ,includeOtherPols);
       		normFactor = GetNormFactorsTM()[l][m][n];
       	}
       	auto it = tField.begin();
       	while (it != tField.end())
       	{
       		if (!isnan(*it))
       		{
       			(*it) *= normFactor;
       		}
       		else
       		{
       			(*it) = 0.;
       		}
       		*it++;
       	}

       	return tField;  // return normalized field.
    }

/*
	double CylindricalCavity::TotalFieldNorm(std::vector<double> field)
	{
		double norm = 0;
		auto it = field.begin();
		while (it != field.end())
		{
			if (!isnan(*it)) norm += (*it)*(*it);
			*it++;
		}
		return sqrt(norm);
	}
*/
	double CylindricalCavity::CalculateDotProductFactor(int bTE, int l, int m, int n, std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, double tThisEventNSamples)
	{
		std::vector<std::vector<std::vector<std::vector<double>>>> tAvgDotProductFactor = GetAvgDotProductFactor();
		tAvgDotProductFactor[bTE][l][m][n] = 1. / ( tThisEventNSamples + 1 ) * ( tAvgDotProductFactor[bTE][l][m][n] * tThisEventNSamples + GetDotProductFactor(tKassParticleXP, anE_normalized, 0) );  // unit velocity \dot unit theta
		SetAvgDotProductFactor(tAvgDotProductFactor);
		return tAvgDotProductFactor[bTE][l][m][n];
	}


    double CylindricalCavity::GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, bool IntermediateFile)
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
    	double unitJdotE = 0.;
    	if ( (tEmag > 0.) && (tVmag > 0.) )
    	{
    		unitJdotE = fabs(tEx*tVx + tEy*tVy)/tEmag/tVmag;
    	}


    	//  Write trajectory points, dot product, and E-field mag to file for debugging etc.
    	if (IntermediateFile)
    	{
            scarab::path dataDir = TOSTRING(PB_DATA_INSTALL_DIR);
            char buffer[60];
            int a = sprintf(buffer, "%s/../output/dotProducts.txt", dataDir.string().c_str());
            const char *fpname = buffer;
            FILE *fp = fopen(fpname, "a");
            fprintf(fp, "%g %g %g %g\n", tKassParticleXP[0], tKassParticleXP[1], unitJdotE, tEmag);
            fclose(fp);

            printf("unitJdotE is %g, r*cos(theta) is %g, r is %g, and theta is %g, eMag is %g\n",
    			unitJdotE, tKassParticleXP[0]*cos(tKassParticleXP[1]), tKassParticleXP[0], tKassParticleXP[1], tEmag); getchar();
    	}

    	return unitJdotE;
    }

    bool CylindricalCavity::InVolume(std::vector<double> tKassParticleXP)
    {
    	double rLocation = tKassParticleXP[0];
    	double zLocation = tKassParticleXP[2];

    	if ((rLocation < GetDimR()) && (fabs(zLocation) < GetDimL()/2.))
    	{
    		return true;
    	}
    	else
    	{
    		return false;
    	}
    }


    void CylindricalCavity::CheckNormalization(int nModes)
    {

       printf("\n \\int{|E_xlm|^2 dV} = \\mu / \\epsilon \\int{|H_xlm|^2 dV} ?\n\n");

    	for (int l=0; l<nModes; l++)
    	{
    		for (int m=1; m<nModes; m++)
    		{
    			for (int n=0; n<nModes; n++)
    			{
    				double normFactor = pow(GetNormFactorsTE()[l][m][n],2.);
    				if (!std::isnan(normFactor)&&(std::isfinite(normFactor)))
    				{
    					printf("TE%d%d%d E %.4g H %.4g\n", l, m, n, Integrate(l,m,n,1,1)*normFactor,
        		    		LMCConst::MuNull()/LMCConst::EpsNull()*Integrate(l,m,n,1,0)*normFactor);
    				}
    				else
    				{
    					printf("TE%d%d%d is undefined.\n", l, m, n);
    				}

    			}
    		}
    	}


    	for (int l=0; l<nModes; l++)
    	{
    		for (int m=1; m<nModes; m++)
    		{
    			for (int n=1; n<nModes; n++)
    			{
    				double normFactor = pow(GetNormFactorsTM()[l][m][n],2.);
    				if (!std::isnan(normFactor)&&(std::isfinite(normFactor)))
    				{
    					printf("TM%d%d%d E %.4g H %.4g\n", l, m, n, Integrate(l,m,n,0,1)*normFactor,
    		    			LMCConst::MuNull()/LMCConst::EpsNull()*Integrate(l,m,n,0,0)*normFactor);
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


    void CylindricalCavity::PrintModeMaps(int nModes, double zSlice, double thetaSlice)
    {
    	scarab::path dataDir = TOSTRING(PB_DATA_INSTALL_DIR);
        std::string sFileName = (dataDir / "../output/ModemapOutput.root").string();

#ifdef ROOT_FOUND

    	FileWriter* aRootHistoWriter = RootHistoWriter::get_instance();
    	aRootHistoWriter->SetFilename(sFileName);
    	aRootHistoWriter->OpenFile("RECREATE");

    	int nbins = GetNPixels();
    	char hbufferEtheta[60]; char hbufferEr[60];
    	char hbufferHtheta[60]; char hbufferHr[60];
    	int a;
    	const char *hname_Etheta = hbufferEtheta;
    	const char *hname_Htheta = hbufferHtheta;
    	const char *hname_Er = hbufferEr;
    	const char *hname_Hr = hbufferHr;

    	char bufferE[60];
    	char bufferH[60];

    	for (int bTE=0; bTE<2; bTE++) // TM and TE
    	{
        	for (int l=0; l<nModes; l++)
        		for (int m=1; m<nModes; m++)
    	    		for (int n=0; n<nModes; n++)
    		    	{
    			    	printf("l m n is %d %d %d\n", l, m, n);
    		    		if (bTE)
    			    	{
    				    	a = sprintf(hbufferEtheta, "TE%d%d%d_Etheta_z%d", l, m, n, (int)(zSlice*1.e3));
    					    a = sprintf(hbufferEr, "TE%d%d%d_Er_z%d", l, m, n, (int)(zSlice*1.e3));
    				    	a = sprintf(hbufferHtheta, "TE%d%d%d_Htheta_z%d", l, m, n, (int)(zSlice*1.e3));
    					    a = sprintf(hbufferHr, "TE%d%d%d_Hr_z%d", l, m, n, (int)(zSlice*1.e3));
    				    }
    				    else
        				{
        					a = sprintf(hbufferEtheta, "TM%d%d%d_Etheta_z%d", l, m, n, (int)(zSlice*1.e3));
    	    				a = sprintf(hbufferEr, "TM%d%d%d_Er_z%d", l, m, n, (int)(zSlice*1.e3));
        					a = sprintf(hbufferHtheta, "TM%d%d%d_Htheta_z%d", l, m, n, (int)(zSlice*1.e3));
    	    				a = sprintf(hbufferHr, "TM%d%d%d_Hr_z%d", l, m, n, (int)(zSlice*1.e3));
    		    		}

    			    	TH2D* hTEtheta = new TH2D(hname_Etheta, hname_Etheta, nbins, -LMCConst::Pi(), LMCConst::Pi(), nbins, 0., GetDimR());
    			 	    TH2D* hTEr = new TH2D(hname_Er, hname_Er, nbins, -LMCConst::Pi(), LMCConst::Pi(), nbins, 0., GetDimR());
    			    	TH2D* hTHtheta = new TH2D(hname_Htheta, hname_Htheta, nbins, -LMCConst::Pi(), LMCConst::Pi(), nbins, 0., GetDimR());
    			 	    TH2D* hTHr = new TH2D(hname_Hr, hname_Hr, nbins, -LMCConst::Pi(), LMCConst::Pi(), nbins, 0., GetDimR());

        				double normFactor = 1.0;
        				if (bTE)
    	    			{
    		    			normFactor = GetNormFactorsTE()[l][m][n];
    			    	}
    	    			else
    		    		{
    			    		normFactor = GetNormFactorsTM()[l][m][n];
        				}
    	    			for (unsigned i=0; i<GetNPixels(); i++)
    		    		{
    			    		double r = ((double)i+0.5)/(GetNPixels())*GetDimR();
    				    	for (unsigned j=0; j<GetNPixels(); j++)
        					{
        						double theta = ((double)j+0.5)/(GetNPixels())*2.*LMCConst::Pi();
            					for (unsigned k=0; k<1; k++)
        	    				{
            	    			    double z = zSlice;
    				    		    std::vector<double> tE;
    					    	    std::vector<double> tH;
        							if (bTE)
    	    						{
    		    						tE = fFieldCore->TE_E(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
    			    					tH = fFieldCore->TE_H(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
    				    			}
        							else
    	    						{
    		    						tE = fFieldCore->TM_E(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
        								tH = fFieldCore->TM_H(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
    	    						}
    		    				    if ((!std::isnan(tE.back())))
    			    			    {
        						        hTEtheta->Fill(theta-LMCConst::Pi(),r,tE.back());
    	    					    }
    		    				    if ((!std::isnan(tE.front())))
    			    			    {
        						    	hTEr->Fill(theta-LMCConst::Pi(),r,tE.front());
        						    }
    		    				    if ((!std::isnan(tH.back())))
    			    			    {
        						        hTHtheta->Fill(theta-LMCConst::Pi(),r,tH.back());
    	    					    }
    		    				    if ((!std::isnan(tH.front())))
    			    			    {
        						    	hTHr->Fill(theta-LMCConst::Pi(),r,tH.front());
        						    }

            					}
        					}
        				}
    	    			aRootHistoWriter->Write2DHisto(hTEtheta);
        				aRootHistoWriter->Write2DHisto(hTEr);
    	    			aRootHistoWriter->Write2DHisto(hTHtheta);
        				aRootHistoWriter->Write2DHisto(hTHr);
        				delete hTEtheta; delete hTEr;
        				delete hTHtheta; delete hTHr;
        			}
    	} // bTE
		aRootHistoWriter->CloseFile();

		PrintModeMapsLongSlice(nModes, 0.);

    	LPROG(lmclog, "\n\nTo plot a mode map:\n"
    			"> root file:output/ModemapOutput.root\n"
    			"# _file0->ls()\n"
    			"# TE011_Etheta->SetLineColor(0)\n"
    			"# TE011_Etheta->SetLineWidth(0)\n"
    			"# TE011_Etheta->DrawCopy(\"pol lego2\")\n"
    			"# TPad *p = (TPad*)c1->cd()\n"
    			"# p->SetTheta(90.); p->SetPhi(0.)\n"
    			"# p->Update()\n"
    			"\n\nMode map files have been generated; press RETURN to continue, or Cntrl-C to quit.");
    	getchar();
#endif

    }

    void CylindricalCavity::PrintModeMapsLongSlice(int nModes, double thetaSlice)
    {
    	scarab::path dataDir = TOSTRING(PB_DATA_INSTALL_DIR);
        std::string sFileName = (dataDir / "../output/ModemapOutput.root").string();

#ifdef ROOT_FOUND

    	FileWriter* aRootHistoWriter = RootHistoWriter::get_instance();
    	aRootHistoWriter->SetFilename(sFileName);
    	aRootHistoWriter->OpenFile("UPDATE");

    	int nbins = GetNPixels();
    	char hbufferEtheta[60]; char hbufferEr[60];
    	char hbufferHtheta[60]; char hbufferHr[60];
    	int a;
    	const char *hname_Etheta = hbufferEtheta;
    	const char *hname_Htheta = hbufferHtheta;
    	const char *hname_Er = hbufferEr;
    	const char *hname_Hr = hbufferHr;

    	char bufferE[60];
    	char bufferH[60];

    	for (int bTE=0; bTE<2; bTE++) // TM and TE
    	{
        	for (int l=0; l<nModes; l++)
        		for (int m=1; m<nModes; m++)
    	    		for (int n=0; n<nModes; n++)
    		    	{
    			    	printf("l m n is %d %d %d\n", l, m, n);
    		    		if (bTE)
    			    	{
    				    	a = sprintf(hbufferEtheta, "TE%d%d%d_Etheta_zLong_theta%d", l, m, n, (int)(thetaSlice*1.e3));
    					    a = sprintf(hbufferEr, "TE%d%d%d_Er_zLong_theta%d", l, m, n, (int)(thetaSlice*1.e3));
    				    	a = sprintf(hbufferHtheta, "TE%d%d%d_Htheta_zLong_theta%d", l, m, n, (int)(thetaSlice*1.e3));
    					    a = sprintf(hbufferHr, "TE%d%d%d_Hr_zLong_theta%d", l, m, n, (int)(thetaSlice*1.e3));
    				    }
    				    else
        				{
        					a = sprintf(hbufferEtheta, "TM%d%d%d_Etheta_thetaLong_r%d", l, m, n, (int)(thetaSlice*1.e3));
    	    				a = sprintf(hbufferEr, "TM%d%d%d_Er_zLong_theta%d", l, m, n, (int)(thetaSlice*1.e3));
        					a = sprintf(hbufferHtheta, "TM%d%d%d_Htheta_zLong_theta%d", l, m, n, (int)(thetaSlice*1.e3));
    	    				a = sprintf(hbufferHr, "TM%d%d%d_Hr_zLong_theta%d", l, m, n, (int)(thetaSlice*1.e3));
    		    		}

    			    	TH2D* hTEtheta = new TH2D(hname_Etheta, hname_Etheta, nbins, -GetDimL()/2., GetDimL()/2., nbins, -GetDimR(), GetDimR());
    			 	    TH2D* hTEr = new TH2D(hname_Er, hname_Er, nbins, -GetDimL()/2., GetDimL()/2., nbins, -GetDimR(), GetDimR());
    			    	TH2D* hTHtheta = new TH2D(hname_Htheta, hname_Htheta, nbins, -GetDimL()/2., GetDimL()/2., nbins, -GetDimR(), GetDimR());
    			 	    TH2D* hTHr = new TH2D(hname_Hr, hname_Hr, nbins, -GetDimL()/2., GetDimL()/2., nbins, -GetDimR(), GetDimR());

        				double normFactor = 1.0;
        				if (bTE)
    	    			{
    		    			normFactor = GetNormFactorsTE()[l][m][n];
    			    	}
    	    			else
    		    		{
    			    		normFactor = GetNormFactorsTM()[l][m][n];
        				}
    	    			for (unsigned i=0; i<GetNPixels(); i++)
    		    		{
    			    		double r = -GetDimR() + ((double)i+0.5)/(GetNPixels())*2.*GetDimR();
    				    	for (unsigned j=0; j<1; j++)
        					{
        						double theta = thetaSlice;
            					for (unsigned k=0; k<GetNPixels(); k++)
        	    				{
            	    			    double z = -GetDimL()/2. + ((double)k+0.5)/(GetNPixels())*GetDimL();
    				    		    std::vector<double> tE;
    					    	    std::vector<double> tH;
        							if (bTE)
    	    						{
    		    						tE = fFieldCore->TE_E(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
    			    					tH = fFieldCore->TE_H(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
    				    			}
        							else
    	    						{
    		    						tE = fFieldCore->TM_E(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
        								tH = fFieldCore->TM_H(GetDimR(),GetDimL(),l,m,n,r,theta,z,0);
    	    						}
    		    				    if ((!std::isnan(tE.back())))
    			    			    {
        						        hTEtheta->Fill(z,r,tE.back());
    	    					    }
    		    				    if ((!std::isnan(tE.front())))
    			    			    {
        						    	hTEr->Fill(z,r,tE.front());
        						    }
    		    				    if ((!std::isnan(tH.back())))
    			    			    {
        						        hTHtheta->Fill(z,r,tH.back());
    	    					    }
    		    				    if ((!std::isnan(tH.front())))
    			    			    {
        						    	hTHr->Fill(z,r,tH.front());
        						    }

            					}
        					}
        				}
    	    			aRootHistoWriter->Write2DHisto(hTEtheta);
        				aRootHistoWriter->Write2DHisto(hTEr);
    	    			aRootHistoWriter->Write2DHisto(hTHtheta);
        				aRootHistoWriter->Write2DHisto(hTHr);
        				delete hTEtheta; delete hTEr;
        				delete hTHtheta; delete hTHr;
        			}
    	} // bTE
		aRootHistoWriter->CloseFile();
#endif

    }


    std::vector<double> CylindricalCavity::GetCavityProbeGain()
    {
    	return fProbeGain;
    }
    void CylindricalCavity::SetCavityProbeGain( double aGain, unsigned index )
    {
    	fProbeGain[index] = aGain;
    }
    std::vector<double> CylindricalCavity::GetCavityProbeZ()
    {
    	return fCavityProbeZ;
    }
    void CylindricalCavity::SetCavityProbeZ ( double aZ, unsigned index )
    {
    	fCavityProbeZ[index] = aZ;
    }
    std::vector<double> CylindricalCavity::GetCavityProbeRFrac()
    {
    	return fCavityProbeRFrac;
    }
    void CylindricalCavity::SetCavityProbeRFrac ( double aFraction, unsigned index )
    {
    	fCavityProbeRFrac[index] = aFraction;
    }
    std::vector<double> CylindricalCavity::GetCavityProbeTheta()
    {
	return fCavityProbeTheta;
    }
    void CylindricalCavity::SetCavityProbeTheta ( double aTheta, unsigned index )
    {
	fCavityProbeTheta[index] = aTheta;
    }


} /* namespace locust */


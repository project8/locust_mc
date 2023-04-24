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
    	fProbeGain( 1.),
		fCavityProbeZ( 0. ),
		fCavityProbeRFrac( 0.5 ),
		fCavityProbePhi( 0.0 ),
		fCavityVolume( 0. )
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

        if ( aParam.has( "cavity-probe-gain" ) )
    	{
    		SetCavityProbeGain(aParam["cavity-probe-gain"]().as_double());
    	}

    	if ( aParam.has( "cavity-probe-z" ) )
    	{
    		SetCavityProbeZ(aParam["cavity-probe-z"]().as_double());
    	}

    	if ( aParam.has( "cavity-probe-r-fraction" ) )
    	{
    		SetCavityProbeRFrac(aParam["cavity-probe-r-fraction"]().as_double());
    	}

     	if ( aParam.has( "cavity-probe-phi" ) )
	{
		SetCavityProbePhi(aParam["cavity-probe-phi"]().as_double());
	}

        /*
                if( aParam.has( "modemap-filename" ) )
                {
                    // TO-DO:  This is where we can plan to read in a mode map.                     *
                    fFieldCore = new TBDModeMapClass();
                }
        */

        fFieldCore = new PozarCylindricalCavity();


        // TO-DO:  Move the next 3 lines to a parent class.
        scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
        fFieldCore->ReadBesselZeroes((dataDir / "BesselZeros.txt").string(), 0 );
        fFieldCore->ReadBesselZeroes((dataDir / "BesselPrimeZeros.txt").string(), 1 );

        SetNormFactorsTE(CalculateNormFactors(GetNModes(),1));
        SetNormFactorsTM(CalculateNormFactors(GetNModes(),0));

        CheckNormalization(GetNModes());  // E fields integrate to 1.0 for both TE and TM modes.
        SetCavityVolume();

        if( aParam.has( "plot-mode-maps" ) )
        {
        	LPROG( lmclog, "If ROOT is available, plotting mode maps to file output/ModeMapOutput.root... " );
        	PrintModeMaps(GetNModes(),1, 0.);
        }

    	return true;
    }

    void CylindricalCavity::SetCavityVolume()
    {
    	fCavityVolume = LMCConst::Pi() * this->GetDimR() * this->GetDimR() * this->GetDimL();
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
            			aModeNormFactor[l][m][n] = 1./Integrate(l,m,n,1,1);
            		}
            		else
            		{
            			aModeNormFactor[l][m][n] = 1./Integrate(l,m,n,0,1);
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

    std::vector<double> CylindricalCavity::GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP)
    {
    	std::vector<double> freqPrime;
    	double vz = tKassParticleXP[5];
    	double term1 = fFieldCore->GetBesselNKPrimeZeros(l,m) / GetDimR();
    	double term2 = n * LMCConst::Pi() / GetDimL();
    	double lambda = 1. / pow( 1. / 4. / LMCConst::Pi() / LMCConst::Pi() * ( term1*term1 + term2*term2 ), 0.5);
    	double lambda_c = 2 * LMCConst::Pi() * GetDimR() / fFieldCore->GetBesselNKPrimeZeros(l,m);
    	double vp = LMCConst::C() / pow( 1. - lambda*lambda/lambda_c/lambda_c, 0.5 );
    	double dopplerShift = vz / vp;
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

    double CylindricalCavity::GetFieldAtProbe(int l, int m, int n, bool includeOtherPols, std::vector<double> tKassParticleXP)
	{
    	double rProbe = this->GetCavityProbeRFrac() * this->GetDimR();
	double phiProbe = this->GetCavityProbePhi();
    	double zProbe = this->GetCavityProbeZ();

	double phiEffective = phiProbe; 
	if(l>0)
	{
		//If mode has phi dependence, mode polarization is set by electron location. Probe coupling must be set relative to that angle
		double phiElectron = tKassParticleXP[1];
		phiEffective = phiProbe - phiElectron;
	}	


		std::vector<double> tProbeLocation = {rProbe, phiEffective, zProbe};
		// Factor of sqrt(fCavityVolume) is being applied to the pre-digitized E-fields
		// to try to reduce dependence of detected power on cavity volume.  This
		// is qualitatively consistent with the volume scaling expected from a cavity
		// experiment, and will also support ongoing normalization studies without
		// necessarily having to retune the digitizer frequently.
		double tEFieldAtProbe = sqrt(fCavityVolume) * GetNormalizedModeField(l,m,n,tProbeLocation,0).back(); //Assumes probe couples to E_theta of mode. If mode is polarized, transforms angle to reference frame of electron
    	return fProbeGain * tEFieldAtProbe;
	}

    std::vector<double> CylindricalCavity::GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP, bool includeOtherPols)
    {
    	double tR = tKassParticleXP[0];
    	double tPhi = tKassParticleXP[1];
    	double tZ = tKassParticleXP[2];
       	std::vector<double> tField;

       	tField = fFieldCore->TE_E(GetDimR(),GetDimL(),l,m,n,tR,tPhi,tZ,1);
       	double normFactor = GetNormFactorsTE()[l][m][n];
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
    	double unitJdotE = fabs(tEx*tVx + tEy*tVy)/tEmag/tVmag;


    	//  Write trajectory points, dot product, and E-field mag to file for debugging etc.
    	if (IntermediateFile)
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

    void CylindricalCavity::CheckNormalization(int nModes)
    {

       printf("\n \\int{|E_xlm|^2 dV} = \\mu / \\epsilon \\int{|H_xlm|^2 dV} ?\n\n");

    	for (int l=0; l<nModes; l++)
    	{
    		for (int m=1; m<nModes; m++)
    		{
    			for (int n=0; n<nModes; n++)
    			{
    				double normFactor = GetNormFactorsTE()[l][m][n] / LMCConst::EpsNull();
    				if (!std::isnan(normFactor)&&(std::isfinite(normFactor)))
    				{
    					printf("TE%d%d%d E %.4g H %.4g\n", l, m, n, LMCConst::EpsNull()*Integrate(l,m,n,1,1)*normFactor,
        		    		LMCConst::MuNull()*Integrate(l,m,n,1,0)*normFactor);
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
    				double normFactor = GetNormFactorsTM()[l][m][n] / LMCConst::EpsNull();
    				if (!std::isnan(normFactor)&&(std::isfinite(normFactor)))
    				{
    					printf("TM%d%d%d E %.4g H %.4g\n", l, m, n, LMCConst::EpsNull()*Integrate(l,m,n,0,1)*normFactor,
    		    			LMCConst::MuNull()*Integrate(l,m,n,0,0)*normFactor);
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



    void CylindricalCavity::PrintModeMaps(int nModes, bool bTE, double zSlice)
    {

#ifdef ROOT_FOUND

	    FileWriter* aRootHistoWriter = RootHistoWriter::get_instance();
	    aRootHistoWriter->SetFilename("output/ModeMapOutput.root");
	    aRootHistoWriter->OpenFile("RECREATE");

    	int nbins = this->GetNPixels();
    	char hbuffertheta[60]; char hbufferr[60]; int a;
    	const char *hname_theta = hbuffertheta;
    	const char *hname_r = hbufferr;

    	char bufferE[60];
    	char bufferH[60];

    	for (int l=0; l<nModes; l++)
    		for (int m=1; m<nModes; m++)
    			for (int n=0; n<nModes; n++)
    			{
    				printf("l m n is %d %d %d\n", l, m, n);
    		    	a = sprintf(hbuffertheta, "TE%d%d%d_Etheta", l, m, n);
    		    	a = sprintf(hbufferr, "TE%d%d%d_Er", l, m, n);
    				TH2D* hTEtheta = new TH2D(hname_theta, hname_theta, nbins, -LMCConst::Pi(), LMCConst::Pi(), nbins, 0., this->GetDimR());
    				TH2D* hTEr = new TH2D(hname_r, hname_r, nbins, -LMCConst::Pi(), LMCConst::Pi(), nbins, 0., this->GetDimR());

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

        					}
    					}
    				}
    				aRootHistoWriter->Write2DHisto(hTEtheta);
    				aRootHistoWriter->Write2DHisto(hTEr);
    				delete hTEtheta; delete hTEr;
    			}
		aRootHistoWriter->CloseFile();
    	LPROG(lmclog, "\n\nTo plot a mode map:\n"
    			"> root file:output/ModeMapOutput.root\n"
    			"# _file0->ls()\n"
    			"# hTEtheta->SetLineColor(0)\n"
    			"# hTEtheta->SetLineWidth(0)\n"
    			"# TE011_Etheta->DrawCopy(\"pol lego2\")\n"
    			"# TPad *p = (TPad*)c1->cd()\n"
    			"# p->SetTheta(90.); p->SetPhi(0.)\n"
    			"# p->Update()\n"
    			"\n\nMode map files have been generated; press RETURN to continue, or Cntrl-C to quit.");
    	getchar();
#endif

    }

    double CylindricalCavity::GetCavityProbeGain()
    {
    	return fProbeGain;
    }
    void CylindricalCavity::SetCavityProbeGain( double aGain )
    {
    	fProbeGain = aGain;
    }
    double CylindricalCavity::GetCavityProbeZ()
    {
    	return fCavityProbeZ;
    }
    void CylindricalCavity::SetCavityProbeZ ( double aZ )
    {
    	fCavityProbeZ = aZ;
    }
    double CylindricalCavity::GetCavityProbeRFrac()
    {
    	return fCavityProbeRFrac;
    }
    void CylindricalCavity::SetCavityProbeRFrac ( double aFraction )
    {
    	fCavityProbeRFrac = aFraction;
    }
    double CylindricalCavity::GetCavityProbePhi()
    {
	return fCavityProbePhi;
    }
    void CylindricalCavity::SetCavityProbePhi ( double aPhi )
    {
	fCavityProbePhi = aPhi;
    }


} /* namespace locust */


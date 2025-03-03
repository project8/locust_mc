/*
 * LMCRectangularWaveguide.cc
 *
 *  Created on: Jun 9, 2021
 *      Author: pslocum
 */

#include "LMCRectangularWaveguide.hh"


namespace locust
{
    LOGGER( lmclog, "RectangularWaveguide" );
    RectangularWaveguide::RectangularWaveguide():
    fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }
    RectangularWaveguide::~RectangularWaveguide() {}


    bool RectangularWaveguide::Configure( const scarab::param_node& aParam)
    {

    	if( !Field::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring Field class from CylindricalCavity subclass");
    		return false;
    	}

    	if( aParam.has( "waveguide-central-frequency" ) )
    	{
    		SetCentralFrequency( 2.*LMCConst::Pi()*aParam["waveguide-central-frequency"]().as_double() );
    	}

        if( aParam.has( "waveguide-x" ) )
        {
            SetDimX( aParam["waveguide-x"]().as_double() );
        }

        if( aParam.has( "waveguide-y" ) )
        {
        	SetDimY( aParam["waveguide-y"]().as_double() );
        }

        if( aParam.has( "waveguide-z") )
        {
        	SetDimL( aParam["waveguide-z"]().as_double() );
        }

    	if( aParam.has( "center-to-short" ) ) // for use in waveguide
        {
            SetCenterToShort(aParam["center-to-short"]().as_double());
        }

        if( aParam.has( "center-to-antenna" ) ) // for use in waveguide
        {
            SetCenterToAntenna(aParam["center-to-antenna"]().as_double());
        }

        fFieldCore = new PozarRectangularWaveguide();

        SetNormFactors(CalculateNormFactors(GetNModes(), 1));
        CheckNormalization(GetNModes(), 1);  // E fields integrate to 1.0 for both TE and TM modes.

        if( PlotModeMaps() )
        {
        	LPROG( lmclog, "If ROOT is available, plotting mode maps to file output/ModemapOutput.root... " );
        	PrintModeMaps(GetNModes(), 0., 0.);
        }

        return true;

    }


    double RectangularWaveguide::Integrate(int l, int m, int n, bool teMode, bool eField)
    {

    	std::vector<double> aField;
    	double xPozar, yPozar, xKass, yKass = 0.;
    	double dX = GetDimX()/GetNPixels();
    	double dY = GetDimY()/GetNPixels();
    	double tArea = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<GetNPixels(); i++)
    	{
    		for (unsigned j=0; j<GetNPixels(); j++)
    		{
    	    	xPozar = (double)i*dX;
    	    	yPozar = (double)j*dY;
   	    		xKass = xPozar - GetDimX()/2.;
   	    		yKass = yPozar - GetDimY()/2.;

    	    	if (teMode)
    	    	{
    	    		if (eField)
    	    		{
    	    		    aField = fFieldCore->TE_E(GetDimX(), GetDimY(), m, n, xKass, yKass, GetCentralFrequency());
    	    		}
    	    		else
    	    		{
    	    			aField = fFieldCore->TE_H(GetDimX(), GetDimY(), m, n, xKass, yKass, GetCentralFrequency());
    	    		}
    	    	}
    	    	else
    	    	{
    	    		if (eField)
    	    		{
    	    			aField = fFieldCore->TM_E(GetDimX(), GetDimY(), m, n, xKass, yKass, GetCentralFrequency());
    	    		}
    	    		else
    	    		{
    	    			aField = fFieldCore->TM_H(GetDimX(), GetDimY(), m, n, xKass, yKass, GetCentralFrequency());
    	    		}
    	    	}

    	    	double aFieldMagSq = 0.;
    	    	auto it = aField.begin();
    	    	while (it != aField.end())
    	    	{
		    		if (!std::isnan(*it))
		    			aFieldMagSq += (*it)*(*it);
    	    		*it++;
    	    	}

    			tIntegral += aFieldMagSq*dX*dY;
    		    tArea += dX*dY;  // sanity check area integral.
    		}
    	}
//    	printf("tArea is %g\n", tArea); getchar();
    	return tIntegral;
    }


    double RectangularWaveguide::GetGroupVelocity(int m, int n, double fcyc)
    {
    	double CutOffFrequency = 0.;
    	if ((m<2)&&(n<1))  // most likely case
    	{
    		// rad/s
    		CutOffFrequency = LMCConst::C() * LMCConst::Pi() / GetDimX();
    	}
    	else  // general case
    	{
    		// rad/s
    		CutOffFrequency = LMCConst::C() *
    				sqrt(pow(m*LMCConst::Pi()/GetDimX(),2.) + sqrt(pow(n*LMCConst::Pi()/GetDimY(),2.)));
    	}
        double GroupVelocity = LMCConst::C() * pow( 1. - pow(CutOffFrequency/fcyc, 2.) , 0.5);
        //        printf("GroupVelocity is %g\n", GroupVelocity); getchar();
        return GroupVelocity;
    }


    std::vector<double> RectangularWaveguide::GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP)
    {
    	std::vector<double> freqPrime;
    	double fcyc = tKassParticleXP[7];
    	double groupVelocity = GetGroupVelocity(m,n,fcyc);
    	double zVelocity = 0.;

    	for (unsigned towardAntenna=0; towardAntenna<2; towardAntenna++)
    	{
			zVelocity = ( 1. - towardAntenna*2. ) * tKassParticleXP[5];
            double gammaZ = 1.0 / pow(1.0-pow(zVelocity/groupVelocity,2.),0.5);
            freqPrime.push_back(fcyc * gammaZ * (1.+zVelocity/groupVelocity) );
    	}
    	return freqPrime;
    }

    double RectangularWaveguide::ScaleEPoyntingVector(double fcyc)
    {
    	// This function calculates the coefficients of the Poynting vector integral
    	// in the TE10 mode in WR42.  It then returns the sqrt of the half of the propagating
    	// power that is moving toward the antenna.
    	// After Pozar p. 114:
    	double k = fcyc / LMCConst::C();
    	double k1 = LMCConst::Pi() / GetDimX();
    	double beta = sqrt( k*k - k1*k1 );
    	double areaIntegral = fcyc * LMCConst::MuNull() * pow(GetDimX(),3.) * GetDimY() * beta / 4. / LMCConst::Pi() / LMCConst::Pi();
    	// Modify Pozar integral to allow for power that flows in both directions:
    	areaIntegral /= 2.;
    	// sqrt of propagating power gives amplitude of E
    	return sqrt(areaIntegral);
    }



    double RectangularWaveguide::Z_TE(int l, int m, int n, double fcyc) const
    {
    	double k1 = m * LMCConst::Pi() / GetDimX();
    	double k2 = n * LMCConst::Pi() / GetDimY();
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );
    	double k = fcyc / LMCConst::C();
    	double beta = sqrt(k*k - kc*kc);

    	double Z_TE = k*eta/beta;  // This is 448 ohms for TE10 at 25.9 GHz.
    	return Z_TE;
    }

    double RectangularWaveguide::Z_TM(int l, int m, int n, double fcyc) const
    {
    	double k1 = m * LMCConst::Pi() / GetDimX();
    	double k2 = n * LMCConst::Pi() / GetDimY();
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );
    	double k = fcyc / LMCConst::C();
    	double beta = sqrt(k*k - kc*kc);

    	double Z_TM = beta*eta/k;
    	return Z_TM;
    }


    std::vector<double> RectangularWaveguide::GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP,  bool includeOtherPols, bool bTE)
    {
     	// The l index is inert in the waveguide.
      	double tX = tKassParticleXP[0] * cos(tKassParticleXP[1]);
      	double tY = tKassParticleXP[0] * sin(tKassParticleXP[1]);
      	double fcyc = tKassParticleXP[7];
        std::vector<double> tTE_E_electron = fFieldCore->TE_E(GetDimX(),GetDimY(),m,n,tX,tY,fcyc);
  		double normFactor = GetNormFactors()[bTE][l][m][n];

  		auto it = tTE_E_electron.begin();
  		while (it != tTE_E_electron.end())
  		{
  			if (!std::isnan(*it))
  			{
  				(*it) *= normFactor;
  			}
  			else
  			{
  				(*it) = 0.;
  			}
  			*it++;
  		}

      	return tTE_E_electron;  // return normalized field.
    }

	double RectangularWaveguide::CalculateDotProductFactor(int l, int m, int n, std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, double tThisEventNSamples)
	{
		std::vector<std::vector<std::vector<double>>> tAvgDotProductFactor = GetAvgDotProductFactor();
		tAvgDotProductFactor[l][m][n] = 1. / ( tThisEventNSamples + 1 ) * ( tAvgDotProductFactor[l][m][n] * tThisEventNSamples + GetDotProductFactor(tKassParticleXP, anE_normalized, 0) );  // unit velocity \dot unit theta
		SetAvgDotProductFactor(tAvgDotProductFactor);
		return tAvgDotProductFactor[l][m][n];
	}


    double RectangularWaveguide::GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized, bool IntermediateFile)
    {
    	double tThetaParticle = tKassParticleXP[1];
    	double tEy = aTE_E_normalized.back();
     	double tEmag = fabs(tEy);
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double tVmag = pow(tVx*tVx + tVy*tVy, 0.5);
    	double unitJdotE = 0.;
    	if ( (tEmag > 0.) && (tVmag > 0.) )
    	{
    		unitJdotE = fabs(0. + tEy*tVy)/tEmag/tVmag;
    	}

    	//  Write trajectory points, dot product, and E-field mag to file for debugging etc.
    	if (IntermediateFile)
    	{
            char cBufferFileName[60];
            int a = sprintf(cBufferFileName, "%s/dotProducts.txt", GetOutputPath().c_str());
            const char *fpname = cBufferFileName;
            FILE *fp = fopen(fpname, "a");
            fprintf(fp, "%g %g %g %g\n", tKassParticleXP[0], tKassParticleXP[1], unitJdotE, tEmag);
            fclose(fp);

            printf("unitJdotE is %g, r*cos(theta) is %g, r is %g, and theta is %g, eMag is %g\n",
    			unitJdotE, tKassParticleXP[0]*cos(tKassParticleXP[1]), tKassParticleXP[0], tKassParticleXP[1], tEmag); getchar();
    	}

    	return unitJdotE;
    }



    bool RectangularWaveguide::InVolume(std::vector<double> tKassParticleXP)
    {
    	// TO-DO:  Check x and y axes of Kass vs. waveguide geometry.  For
    	// now, rely on zLocation to define waveguide volume, and assume
    	// that the waveguide jacket termination works to stop particles
    	// when they hit the walls of the waveguide.

//    	double xLocation = tKassParticleXP[0] * cos(tKassParticleXP[1]);
//    	double yLocation = tKassParticleXP[0] * sin(tKassParticleXP[1]);
    	double zLocation = tKassParticleXP[2];

    	if (fabs(zLocation) < GetDimL()/2.)
    	{
    		return true;
    	}
    	else
    	{
    		return false;
    	}
    }


    void RectangularWaveguide::PrintModeMaps(int nModes, double zSlice, double thetaSlice)
    {
        std::string sFileName = GetOutputPath()+"/ModemapOutput.root";
        bool bTE = 1; // TE and not TM.

#ifdef ROOT_FOUND
        FileWriter* aRootHistoWriter = RootHistoWriter::get_instance();
        aRootHistoWriter->SetFilename(sFileName);
        aRootHistoWriter->SetFilename(sFileName);
        aRootHistoWriter->OpenFile("RECREATE");

        int nbins = GetNPixels();
        char hbufferx[60]; char hbuffery[60]; int a;
        char labelx[60]; char labely[60];

        for (int l=0; l<1; l++)
        	for (int m=1; m<nModes; m++)
        		for (int n=0; n<nModes; n++)
        		{
        			printf("l m n is %d %d %d\n", l, m, n);
        			a = sprintf(hbuffery, "TE%d%d_Ey", m, n);
        			const char *hname_y = hbuffery;
        			a = sprintf(labely, "TE%d%d_Ey; x (m); y (m)", m, n);
        			const char *haxis_y = labely;
        			a = sprintf(hbufferx, "TE%d%d_Ex", m, n);
        			const char *hname_x = hbufferx;
        			a = sprintf(labelx, "TE%d%d_Ex; x (m); y (m)", m, n);
        			const char *haxis_x = labelx;
        			TH2D* hTEy = new TH2D(hname_y, haxis_y, nbins, -GetDimX()/2., GetDimX()/2., nbins, -GetDimY()/2., GetDimY()/2.);
        			TH2D* hTEx = new TH2D(hname_x, haxis_x, nbins, -GetDimX()/2., GetDimX()/2., nbins, -GetDimY()/2., GetDimY()/2.);

        			double normFactor = GetNormFactors()[bTE][l][m][n];
        			int a = 0;

        			for (unsigned i=0; i<GetNPixels(); i++)
        			{
        				double x = ((double)i+0.5)/GetNPixels()*GetDimX() - GetDimX()/2.;
        				for (unsigned j=0; j<GetNPixels(); j++)
        				{
        					double y = ((double)j+0.5)/GetNPixels()*GetDimY() - GetDimY()/2.;
        					for (unsigned k=0; k<1; k++)
        					{
        						double z = zSlice;
        						std::vector<double> tE;
        						std::vector<double> tH;
        						tE = fFieldCore->TE_E(GetDimX(),GetDimY(),m,n,x,y,GetCentralFrequency());
        						if (!std::isnan(tE.back())) hTEy->Fill(x,y,tE.back());
        						if (!std::isnan(tE.front())) hTEx->Fill(x,y,tE.front());
        					}
        				}
        			}
        			aRootHistoWriter->Write2DHisto(hTEy);
        			aRootHistoWriter->Write2DHisto(hTEx);
        			delete hTEy; delete hTEx;
        		}

        aRootHistoWriter->CloseFile();
        LPROG(lmclog, "\n\nTo plot a mode map:\n"
    			"> root file:output/ModemapOutput.root\n"
    			"# _file0->ls()\n"
    			"# TE10_Ey->SetLineColor(0)\n"
    			"# TE10_Ey->SetLineWidth(0)\n"
    			"# TE10_Ey->DrawCopy(\"colz\")\n"
    			"\n\nMode map files have been generated; press RETURN to continue, or Cntrl-C to quit.");
        getchar();

#endif

    }



} /* namespace locust */


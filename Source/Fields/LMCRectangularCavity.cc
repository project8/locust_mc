/*
 * LMCRectangularCavity.cc
 *
 *  Created on: Jun 9, 2021
 *      Author: pslocum
 */

#include "LMCRectangularCavity.hh"


namespace locust
{

    LOGGER( lmclog, "RectangularCavity" );
    RectangularCavity::RectangularCavity():
    	fProbeGain( {1., 1.}),
		fCavityProbeZ( {0., 0.} ),
		fCavityProbeRFrac( {0.5, 0.5} ),
		fCavityProbeTheta( {0.0, 0.0} )
		{}

    RectangularCavity::~RectangularCavity() {}



    bool RectangularCavity::Configure( const scarab::param_node& aParam)
    {

    	if( !Field::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring Field class from RectangularCavity subclass");
    		return false;
    	}

        if( aParam.has( "cavity-x" ) )
        {
            SetDimX( aParam["cavity-x"]().as_double() );
        }

        if( aParam.has( "cavity-y" ) )
        {
            SetDimY( aParam["cavity-y"]().as_double() );
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

/*
     	if( aParam.has( "modemap-filename" ) )
     	{
     		fFieldCore = new ModeMapRectangularCavity();
     		if (!fFieldCore->ReadModeMapTE_E(aParam["modemap-filename"]().as_string()))
     		{
     			LERROR(lmclog,"There was a problem uploading the mode map.");
     			exit(-1);
     		}
     	}
     	else
*/
     	{
     		fFieldCore = new PozarRectangularCavity();
     	}

        SetNormFactorsTE(CalculateNormFactors(GetNModes(),1));
        SetNormFactorsTM(CalculateNormFactors(GetNModes(),0));

        CheckNormalization(GetNModes());  // E fields integrate to 1.0 for both TE and TM modes.

        if( PlotModeMaps() )
        {
        	double zSlice = 0.0;
        	if (aParam.has( "map-z-slice" )) zSlice = aParam["map-z-slice"]().as_double();
        	LPROG( lmclog, "If ROOT is available, plotting mode maps to file output/ModemapOutput*.root... " );
        	PrintModeMaps(GetNModes(), zSlice);
        }

    	return true;
    }

    std::vector<std::vector<std::vector<double>>> RectangularCavity::CalculateNormFactors(int nModes, bool bTE)
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



    double RectangularCavity::Integrate(int l, int m, int n, bool teMode, bool eField)
    {

    	std::vector<double> aField;
    	double xPozar, yPozar, zPozar, xKass, yKass, zKass = 0.;
    	double dX = GetDimX()/GetNPixels();
    	double dY = GetDimY()/GetNPixels();
    	double dZ = GetDimL()/GetNPixels();
    	double tVolume = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<GetNPixels(); i++)
    	{
    		for (unsigned j=0; j<GetNPixels(); j++)
    		{
    			for (unsigned k=0; k<GetNPixels(); k++)
    			{

    	    		xPozar = (double)i*dX;
    	    		yPozar = (double)j*dY;
    	    		zPozar = (double)k*dZ;
    	    		xKass = xPozar - GetDimX()/2.;
    	    		yKass = yPozar - GetDimY()/2.;
    	    		zKass = zPozar - GetDimL()/2.;

    	    		if (teMode)
    	    		{
    	    			if (eField)
    	    			{
    	    		    	aField = fFieldCore->TE_E(GetDimX(), GetDimY(), GetDimL(), l, m, n, xKass, yKass, zKass);
    	    			}
    	    			else
    	    			{
    	    				aField = fFieldCore->TE_H(GetDimX(), GetDimY(), GetDimL(), l, m, n, xKass, yKass, zKass);
    	    			}
    	    		}
    	    		else
    	    		{
    	    			if (eField)
    	    			{
    	    				aField = fFieldCore->TM_E(GetDimX(), GetDimY(), GetDimL(), l, m, n, xKass, yKass, zKass);
    	    			}
    	    			else
    	    			{
    	    				aField = fFieldCore->TM_H(GetDimX(), GetDimY(), GetDimL(), l, m, n, xKass, yKass, zKass);
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

    				tIntegral += aFieldMagSq*dX*dY*dZ;
    		    	//tVolume += dX*dY*dZ;  // sanity check volume integral.
    			}
    		}
    	}
    	return tIntegral;
    }

    std::vector<double> RectangularCavity::GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP)
    {
    	// Following the calculations in https://3.basecamp.com/3700981/buckets/3107037/documents/4496563487#__recording_4496575215
    	std::vector<double> freqPrime;

    	double vz = tKassParticleXP[5];
    	double k1 = l * LMCConst::Pi() / GetDimX();
    	double k2 = m * LMCConst::Pi() / GetDimY();
    	double kc = sqrt( k1*k1 + k2*k2 );
    	double lambda_c = 2. * LMCConst::Pi() / kc;
    	double lambda = 1. / sqrt( 1. / lambda_c / lambda_c + n*n / 4. / GetDimL() / GetDimL() );
    	double vp = LMCConst::C() * sqrt( 1. - lambda*lambda / lambda_c/lambda_c );
    	double dopplerShift = 0.;
    	if (vp > 0.) dopplerShift = vz / vp;
		freqPrime.push_back( ( 1. + dopplerShift ) * tKassParticleXP[7] );

    	return freqPrime;
    }


    double RectangularCavity::Z_TE(int l, int m, int n, double fcyc) const
    {

     	double Z_TE = 500.; // ohms
    	return Z_TE;

    }

    double RectangularCavity::Z_TM(int l, int m, int n, double fcyc) const
    {
    	double Z_TM = 500.; //ohms
    	return Z_TM;
    }


    std::vector<double> RectangularCavity::GetTE_E(int l, int m, int n, double x, double y, double z, bool includeOtherPols)
    {
    	return fFieldCore->TE_E(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
    }

    std::vector<double> RectangularCavity::GetTM_E(int l, int m, int n, double x, double y, double z, bool includeOtherPols)
    {
    	return fFieldCore->TM_E(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
    }

    std::vector<double> RectangularCavity::GetFieldAtProbe(int l, int m, int n, bool includeOtherPols, std::vector<double> tKassParticleXP, bool teMode)
    {

    	/* Even for a rectangular cavity, probe location is defined in polar coordinates for consistency
    	 * with LMCKassCurrentTransmitter output. */

    	std::vector<double> rProbe;
        rProbe.push_back(GetCavityProbeRFrac()[0] * GetDimR());
        rProbe.push_back(GetCavityProbeRFrac()[1] * GetDimR());

    	std::vector<double> thetaProbe = GetCavityProbeTheta();
    	std::vector<double> zProbe = GetCavityProbeZ();

    	std::vector<std::vector<double>> tProbeLocation;
    	tProbeLocation.push_back({rProbe[0], thetaProbe[0], zProbe[0]});
    	tProbeLocation.push_back({rProbe[1], thetaProbe[1], zProbe[1]});

    	std::vector<double> tEFieldAtProbe;
    	tEFieldAtProbe.push_back( NormalizedEFieldMag(GetNormalizedModeField(l,m,n,tProbeLocation[0],0,teMode)) );
    	tEFieldAtProbe.push_back( NormalizedEFieldMag(GetNormalizedModeField(l,m,n,tProbeLocation[1],0,teMode)) );

    	return {fProbeGain[0] * tEFieldAtProbe[0], fProbeGain[1] * tEFieldAtProbe[1]};

    }

    std::vector<double> RectangularCavity::GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP, bool includeOtherPols, bool teMode)
    {
    	double tR = tKassParticleXP[0];
    	double tTheta = tKassParticleXP[1];
    	double tZ = tKassParticleXP[2];
    	double tX = tR*cos(tTheta);
    	double tY = tR*sin(tTheta);
       	std::vector<double> tField;
       	double normFactor;
       	if(teMode)
       	{
       		tField = fFieldCore->TE_E(GetDimX(),GetDimY(),GetDimL(),l,m,n,tX,tY,tZ);
       		normFactor = GetNormFactorsTE()[l][m][n];
       	}
       	else
       	{
       		tField = fFieldCore->TM_E(GetDimX(),GetDimY(),GetDimL(),l,m,n,tX,tY,tZ);
       		normFactor = GetNormFactorsTM()[l][m][n];
       	}
       	auto it = tField.begin();
       	while (it != tField.end())
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

       	return tField;  // return normalized field.
    }

	double RectangularCavity::CalculateDotProductFactor(int l, int m, int n, std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, double tThisEventNSamples)
	{
		std::vector<std::vector<std::vector<double>>> tAvgDotProductFactor = GetAvgDotProductFactor();
		tAvgDotProductFactor[l][m][n] = 1. / ( tThisEventNSamples + 1 ) * ( tAvgDotProductFactor[l][m][n] * tThisEventNSamples + GetDotProductFactor(tKassParticleXP, anE_normalized, 0) );  // unit velocity \dot unit theta
		SetAvgDotProductFactor(tAvgDotProductFactor);
		return tAvgDotProductFactor[l][m][n];
	}


    double RectangularCavity::GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, bool IntermediateFile)
    {
    	double tThetaParticle = tKassParticleXP[1];
    	double tEx = 0.;
    	double tEy = 0.;

    	if (!std::isnan(anE_normalized.back()))
    	{
    		tEy = anE_normalized.back();
    	}
    	if (!std::isnan(anE_normalized.front()))
    	{
    		tEx = anE_normalized.front();
    	}
    	double tEmag = pow(tEx*tEx + tEy*tEy, 0.5);
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

    bool RectangularCavity::InVolume(std::vector<double> tKassParticleXP)
    {
    	double xLocation = tKassParticleXP[0]*cos(tKassParticleXP[1]);
    	double yLocation = tKassParticleXP[0]*sin(tKassParticleXP[1]);
    	double zLocation = tKassParticleXP[2];

    	if ((fabs(xLocation) < GetDimX()/2.) && (fabs(yLocation) < GetDimY()/2.) && (fabs(zLocation) < GetDimL()/2.))
    	{
    		return true;
    	}
    	else
    	{
    		return false;
    	}
    }


    void RectangularCavity::CheckNormalization(int nModes)
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


    void RectangularCavity::PrintModeMaps(int nModes, double zSlice)
    {
    	scarab::path dataDir = TOSTRING(PB_DATA_INSTALL_DIR);
        std::string sFileName = (dataDir / "../output/ModemapOutput.root").string();

#ifdef ROOT_FOUND

    	FileWriter* aRootHistoWriter = RootHistoWriter::get_instance();
    	aRootHistoWriter->SetFilename(sFileName);
    	aRootHistoWriter->OpenFile("RECREATE");

    	int nbins = GetNPixels();
    	char hbufferEx[60]; char hbufferEy[60];
    	char hbufferHx[60]; char hbufferHy[60];
    	int a;
    	const char *hname_Ex = hbufferEx;
    	const char *hname_Hx = hbufferHx;
    	const char *hname_Ey = hbufferEy;
    	const char *hname_Hy = hbufferHy;

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
    				    	a = sprintf(hbufferEx, "TE%d%d%d_Ex_z%dmm", l, m, n, (int)(zSlice*1.e3));
    				    	a = sprintf(hbufferEy, "TE%d%d%d_Ey_z%dmm", l, m, n, (int)(zSlice*1.e3));
    				    	a = sprintf(hbufferHx, "TE%d%d%d_Hx_z%dmm", l, m, n, (int)(zSlice*1.e3));
    				    	a = sprintf(hbufferHy, "TE%d%d%d_Hy_z%dmm", l, m, n, (int)(zSlice*1.e3));
    				    }
    				    else
        				{
    				    	a = sprintf(hbufferEx, "TM%d%d%d_Ex_z%dmm", l, m, n, (int)(zSlice*1.e3));
    				    	a = sprintf(hbufferEy, "TM%d%d%d_Ey_z%dmm", l, m, n, (int)(zSlice*1.e3));
    				    	a = sprintf(hbufferHx, "TM%d%d%d_Hx_z%dmm", l, m, n, (int)(zSlice*1.e3));
    				    	a = sprintf(hbufferHy, "TM%d%d%d_Hy_z%dmm", l, m, n, (int)(zSlice*1.e3));
    		    		}

    			    	TH2D* hTEx = new TH2D(hname_Ex, (std::string(hname_Ex)+";x(m);y(m)").c_str(), nbins, -GetDimX()/2., GetDimX()/2., nbins, -GetDimY()/2., GetDimY()/2.);
    			    	TH2D* hTEy = new TH2D(hname_Ey, (std::string(hname_Ey)+";x(m);y(m)").c_str(), nbins, -GetDimX()/2., GetDimX()/2., nbins, -GetDimY()/2., GetDimY()/2.);
    			    	TH2D* hTHx = new TH2D(hname_Hx, (std::string(hname_Hx)+";x(m);y(m)").c_str(), nbins, -GetDimX()/2., GetDimX()/2., nbins, -GetDimY()/2., GetDimY()/2.);
    			    	TH2D* hTHy = new TH2D(hname_Hy, (std::string(hname_Hy)+";x(m);y(m)").c_str(), nbins, -GetDimX()/2., GetDimX()/2., nbins, -GetDimY()/2., GetDimY()/2.);

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
    			    		double x = ((double)i+0.5)/(GetNPixels())*GetDimX() - GetDimX()/2.;
    				    	for (unsigned j=0; j<GetNPixels(); j++)
        					{
        						double y = ((double)j+0.5)/(GetNPixels())*GetDimY() - GetDimY()/2.;
            					for (unsigned k=0; k<1; k++)
        	    				{
            	    			    double z = zSlice;
    				    		    std::vector<double> tE;
    					    	    std::vector<double> tH;
        							if (bTE)
    	    						{
    		    						tE = fFieldCore->TE_E(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
    			    					tH = fFieldCore->TE_H(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
    				    			}
        							else
    	    						{
    		    						tE = fFieldCore->TM_E(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
        								tH = fFieldCore->TM_H(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
    	    						}
    		    				    if ((!std::isnan(tE.back())))
    			    			    {
        						        hTEy->Fill(x,y,tE.back());
    	    					    }
    		    				    if ((!std::isnan(tE.front())))
    			    			    {
        						    	hTEx->Fill(x,y,tE.front());
        						    }
    		    				    if ((!std::isnan(tH.back())))
    			    			    {
        						        hTHy->Fill(x,y,tH.back());
    	    					    }
    		    				    if ((!std::isnan(tH.front())))
    			    			    {
        						    	hTHx->Fill(x,y,tH.front());
        						    }

            					}
        					}
        				}
    	    			aRootHistoWriter->Write2DHisto(hTEx);
        				aRootHistoWriter->Write2DHisto(hTEy);
    	    			aRootHistoWriter->Write2DHisto(hTHx);
        				aRootHistoWriter->Write2DHisto(hTHy);
        				delete hTEx; delete hTEy;
        				delete hTHx; delete hTHy;
        			}
    	} // bTE
		aRootHistoWriter->CloseFile();

		PrintModeMapsLongSlice(nModes, GetDimX()*0.1);

    	LPROG(lmclog, "\n\nTo plot a mode map:\n"
    			"> root file:output/ModemapOutput.root\n"
    			"# _file0->ls()\n"
    			"# TE011->DrawCopy()\n"
    			"\n\nMode map files have been generated; press RETURN to continue, or Cntrl-C to quit.");
    	getchar();
#endif

    }

    void RectangularCavity::PrintModeMapsLongSlice(int nModes, double xSlice)
    {
    	scarab::path dataDir = TOSTRING(PB_DATA_INSTALL_DIR);
        std::string sFileName = (dataDir / "../output/ModemapOutput.root").string();

#ifdef ROOT_FOUND

    	FileWriter* aRootHistoWriter = RootHistoWriter::get_instance();
    	aRootHistoWriter->SetFilename(sFileName);
    	aRootHistoWriter->OpenFile("UPDATE");

    	int nbins = GetNPixels();
    	char hbufferEx[60]; char hbufferEy[60];
    	char hbufferHx[60]; char hbufferHy[60];
    	int a;
    	const char *hname_Ex = hbufferEx;
    	const char *hname_Hx = hbufferHx;
    	const char *hname_Ey = hbufferEy;
    	const char *hname_Hy = hbufferHy;

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
    				    	a = sprintf(hbufferEx, "TE%d%d%d_Ex_zLong_x%dmm", l, m, n, (int)(xSlice*1000.));
    				    	a = sprintf(hbufferEy, "TE%d%d%d_Ey_zLong_x%dmm", l, m, n, (int)(xSlice*1000.));
    				    	a = sprintf(hbufferHx, "TE%d%d%d_Hx_zLong_x%dmm", l, m, n, (int)(xSlice*1000.));
    				    	a = sprintf(hbufferHy, "TE%d%d%d_Hy_zLong_x%dmm", l, m, n, (int)(xSlice*1000.));
    				    }
    				    else
        				{
    				    	a = sprintf(hbufferEx, "TM%d%d%d_Ex_zLong_x%dmm", l, m, n, (int)(xSlice*1000.));
    				    	a = sprintf(hbufferEy, "TM%d%d%d_Ey_zLong_x%dmm", l, m, n, (int)(xSlice*1000.));
    				    	a = sprintf(hbufferHx, "TM%d%d%d_Hx_zLong_x%dmm", l, m, n, (int)(xSlice*1000.));
    				    	a = sprintf(hbufferHy, "TM%d%d%d_Hy_zLong_x%dmm", l, m, n, (int)(xSlice*1000.));
    		    		}

    		    		TH2D* hTEx = new TH2D(hname_Ex, (std::string(hname_Ex)+";z(m);y(m)").c_str(), nbins, -GetDimL()/2., GetDimL()/2., nbins, -GetDimY()/2., GetDimY()/2.);
    		    		TH2D* hTEy = new TH2D(hname_Ey, (std::string(hname_Ey)+";z(m);y(m)").c_str(), nbins, -GetDimL()/2., GetDimL()/2., nbins, -GetDimY()/2., GetDimY()/2.);
    		    		TH2D* hTHx = new TH2D(hname_Hx, (std::string(hname_Hx)+";z(m);y(m)").c_str(), nbins, -GetDimL()/2., GetDimL()/2., nbins, -GetDimY()/2., GetDimY()/2.);
    		    		TH2D* hTHy = new TH2D(hname_Hy, (std::string(hname_Hy)+";z(m);y(m)").c_str(), nbins, -GetDimL()/2., GetDimL()/2., nbins, -GetDimY()/2., GetDimY()/2.);

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
    			    		double x = xSlice;
    				    	for (unsigned j=0; j<1; j++)
        					{
        			    		double y = -GetDimY()/2. + ((double)i+0.5)/(GetNPixels())*GetDimY();
            					for (unsigned k=0; k<GetNPixels(); k++)
        	    				{
            	    			    double z = -GetDimL()/2. + ((double)k+0.5)/(GetNPixels())*GetDimL();
    				    		    std::vector<double> tE;
    					    	    std::vector<double> tH;
        							if (bTE)
    	    						{
    		    						tE = fFieldCore->TE_E(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
    			    					tH = fFieldCore->TE_H(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
    				    			}
        							else
    	    						{
    		    						tE = fFieldCore->TM_E(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
        								tH = fFieldCore->TM_H(GetDimX(),GetDimY(),GetDimL(),l,m,n,x,y,z);
    	    						}
    		    				    if ((!std::isnan(tE.back())))
    			    			    {
        						        hTEy->Fill(z,y,tE.back());
    	    					    }
    		    				    if ((!std::isnan(tE.front())))
    			    			    {
        						    	hTEx->Fill(z,y,tE.front());
        						    }
    		    				    if ((!std::isnan(tH.back())))
    			    			    {
        						        hTHx->Fill(z,y,tH.back());
    	    					    }
    		    				    if ((!std::isnan(tH.front())))
    			    			    {
        						    	hTHx->Fill(z,y,tH.front());
        						    }

            					}
        					}
        				}
    	    			aRootHistoWriter->Write2DHisto(hTEx);
        				aRootHistoWriter->Write2DHisto(hTEy);
    	    			aRootHistoWriter->Write2DHisto(hTHx);
        				aRootHistoWriter->Write2DHisto(hTHy);
        				delete hTEx; delete hTEy;
        				delete hTHx; delete hTHy;
        			}
    	} // bTE
		aRootHistoWriter->CloseFile();
#endif

    }


    std::vector<double> RectangularCavity::GetCavityProbeGain()
    {
    	return fProbeGain;
    }
    void RectangularCavity::SetCavityProbeGain( double aGain, unsigned index )
    {
    	fProbeGain[index] = aGain;
    }
    std::vector<double> RectangularCavity::GetCavityProbeZ()
    {
    	return fCavityProbeZ;
    }
    void RectangularCavity::SetCavityProbeZ ( double aZ, unsigned index )
    {
    	fCavityProbeZ[index] = aZ;
    }
    std::vector<double> RectangularCavity::GetCavityProbeRFrac()
    {
    	return fCavityProbeRFrac;
    }
    void RectangularCavity::SetCavityProbeRFrac ( double aFraction, unsigned index )
    {
    	fCavityProbeRFrac[index] = aFraction;
    }
    std::vector<double> RectangularCavity::GetCavityProbeTheta()
    {
	return fCavityProbeTheta;
    }
    void RectangularCavity::SetCavityProbeTheta ( double aTheta, unsigned index )
    {
	fCavityProbeTheta[index] = aTheta;
    }


} /* namespace locust */


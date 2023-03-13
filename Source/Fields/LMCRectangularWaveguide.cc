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

    	LWARN( lmclog, "The rectanguar waveguide simulation is not fully supported with a suitable trap right now ... " );
    	LWARN( lmclog, "Press return to continue ... " );
    	getchar();

    	if( !Field::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring Field class from CylindricalCavity subclass");
    		return false;
    	}

    	if( aParam.has( "center-to-short" ) ) // for use in e-gun
        {
            fInterface->fCENTER_TO_SHORT = aParam["center-to-short"]().as_double();
        }

        if( aParam.has( "center-to-antenna" ) ) // for use in e-gun
        {
            fInterface->fCENTER_TO_ANTENNA = aParam["center-to-antenna"]().as_double();
        }

        fFieldCore = new PozarRectangular();

        SetNormFactorsTE(CalculateNormFactors(GetNModes(),1));
        SetNormFactorsTM(CalculateNormFactors(GetNModes(),0));

        CheckNormalization(GetNModes());  // E fields integrate to 1.0 for both TE and TM modes.

        if( aParam.has( "mode-maps" ) )
        {
        	PrintModeMaps(GetNModes(),1);
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


    std::vector<double> PozarRectangular::TE_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc)
    {

    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;

    	// from Pozar
    	std::vector<double> TE_E;
    	double k1 = m * LMCConst::Pi() / dimX;
    	double k2 = n * LMCConst::Pi() / dimY;
    	double kc = pow(k1*k1+k2*k2,0.5);

    	double tEx = fcyc*LMCConst::MuNull()*n*LMCConst::Pi()/kc/kc/dimY * cos(k1*x) * sin(k2*y);
    	double tEy = -fcyc*LMCConst::MuNull()*m*LMCConst::Pi()/kc/kc/dimX * sin(k1*x) * cos(k2*y);

    	TE_E.push_back(tEx);
    	TE_E.push_back(tEy);
        return TE_E;
    }


    std::vector<double> PozarRectangular::TE_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc)
    {
    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;

    	// from Pozar
    	std::vector<double> TE_H;
    	double k1 = m * LMCConst::Pi() / dimX;
    	double k2 = n * LMCConst::Pi() / dimY;
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double k = fcyc * sqrt(LMCConst::EpsNull()*LMCConst::MuNull());
    	double beta = sqrt(k*k - kc*kc);

    	double tHx = beta*m*LMCConst::Pi()/kc/kc/dimX * sin(k1*x) * cos(k2*y);
    	double tHy = beta*n*LMCConst::Pi()/kc/kc/dimY * cos(k1*x) * sin(k2*y);

    	TE_H.push_back(tHx);
    	TE_H.push_back(tHy);
        return TE_H;
    }


    std::vector<double> PozarRectangular::TM_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc)
    {
    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;

    	// from Pozar
    	std::vector<double> TM_E;
    	double k1 = m * LMCConst::Pi() / dimX;
    	double k2 = n * LMCConst::Pi() / dimY;
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double k = fcyc * sqrt(LMCConst::EpsNull()*LMCConst::MuNull());
    	double beta = sqrt(k*k - kc*kc);

    	double tEx = beta*m*LMCConst::Pi()/kc/kc/dimX * cos(k1*x) * sin(k2*y);
    	double tEy = beta*n*LMCConst::Pi()/kc/kc/dimY * sin(k1*x) * cos(k2*y);
    	TM_E.push_back(tEx);
    	TM_E.push_back(tEy);
        return TM_E;
    }

    std::vector<double> PozarRectangular::TM_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc)
    {
    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;

    	// from Pozar
    	std::vector<double> TM_H;
    	double k1 = m * LMCConst::Pi() / dimX;
    	double k2 = n * LMCConst::Pi() / dimY;
    	double kc = pow(k1*k1+k2*k2,0.5);

    	double tHx = fcyc*LMCConst::EpsNull()*n*LMCConst::Pi()/kc/kc/dimY * sin(k1*x) * cos(k2*y);
    	double tHy = fcyc*LMCConst::EpsNull()*m*LMCConst::Pi()/kc/kc/dimX * cos(k1*x) * sin(k2*y);
    	TM_H.push_back(tHx);
    	TM_H.push_back(tHy);
        return TM_H;
    }

    std::vector<double> RectangularWaveguide::GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP)
    {
     	// The l index is inert in the waveguide.
      	double tX = tKassParticleXP[0] * cos(tKassParticleXP[1]);
      	double tY = tKassParticleXP[0] * sin(tKassParticleXP[1]);
      	double fcyc = tKassParticleXP[7];
        std::vector<double> tTE_E_electron = fFieldCore->TE_E(GetDimX(),GetDimY(),m,n,tX,tY,fcyc);
  		double normFactor = GetNormFactorsTE()[l][m][n];

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


    double RectangularWaveguide::GetDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized, bool IntermediateFile)
    {
    	double tThetaParticle = tKassParticleXP[1];
    	double tEy = aTE_E_normalized.back();
     	double tEmag = fabs(tEy);
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double tVmag = pow(tVx*tVx + tVy*tVy, 0.5);
    	double unitJdotE = fabs(0. + tEy*tVy)/tEmag/tVmag;


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

    std::vector<std::vector<std::vector<double>>> RectangularWaveguide::CalculateNormFactors(int nModes, bool bTE)
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


    void RectangularWaveguide::CheckNormalization(int nModes)
    {

        printf("\n |E_mn|^2 dA = 1.0.  |H_mn| can vary.  Index l is not used in the waveguide.\n");
        printf("m is the index in the x-direction (widest).  n is the index in the y-direction (narrowest).\n");
        printf("The waveguide calculations assume a signal frequency of 25.9e9 Hz. This can be changed "
        		"by adjusting the parameter \"central-frequency\" on the command line.\n\n");


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



    void RectangularWaveguide::PrintModeMaps(int nModes, bool bTE)
    {

    	char bufferE[60];
    	char bufferH[60];
    	unsigned modeCounter = 0;

    	for (int l=0; l<nModes; l++)
    		for (int m=1; m<nModes; m++)
    			for (int n=0; n<nModes; n++)
    			{
    				printf("l m n is %d %d %d\n", l, m, n);
    				double normFactor = 1.0;
    				int a = 0;
    				if (bTE)
    				{
    					normFactor = GetNormFactorsTE()[l][m][n];
        				a = sprintf(bufferE, "output/ModeMapTE%d%d%d_E.txt", l, m, n);
        				a = sprintf(bufferH, "output/ModeMapTE%d%d%d_H.txt", l, m, n);
    				}
    				else
    				{
    					normFactor = GetNormFactorsTM()[l][m][n];
        				a = sprintf(bufferE, "output/ModeMapTM%d%d%d_E.txt", l, m, n);
        				a = sprintf(bufferH, "output/ModeMapTM%d%d%d_H.txt", l, m, n);
    				}
    				const char *fpnameE = bufferE;
    				FILE *fp_E = fopen(fpnameE, "w");
    				const char *fpnameH = bufferH;
    				FILE *fp_H = fopen(fpnameH, "w");
    				for (unsigned i=0; i<GetNPixels()+1; i++)
    				{
    					double x = (double)i/GetNPixels()*GetDimX() - GetDimX()/2.;
    					for (unsigned j=0; j<GetNPixels()+1; j++)
    					{
        					double y = (double)j/GetNPixels()*GetDimY() - GetDimY()/2.;
        					for (unsigned k=0; k<GetNPixels()+1; k++)
        					{
            				    double z = (double)k/GetNPixels()*GetDimL() - GetDimL()/2.;
    						    std::vector<double> tE;
    						    std::vector<double> tH;
    							tE = fFieldCore->TE_E(GetDimX(),GetDimY(),m,n,x,y,GetCentralFrequency());
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



} /* namespace locust */


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
    CylindricalCavity::CylindricalCavity() {}

    CylindricalCavity::~CylindricalCavity() {}

    PozarCylindrical::PozarCylindrical() {}
    PozarCylindrical::~PozarCylindrical() {}



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

        fFieldCore = new PozarCylindrical();
        scarab::path dataDir = aParam.get_value( "data-dir", ( TOSTRING(PB_DATA_INSTALL_DIR) ) );
        fFieldCore->ReadBesselZeroes((dataDir / "BesselZeros.txt").string(), 0 );
        fFieldCore->ReadBesselZeroes((dataDir / "BesselPrimeZeros.txt").string(), 1 );

        SetNormFactorsTE(CalculateNormFactors(GetNModes(),1));
        SetNormFactorsTM(CalculateNormFactors(GetNModes(),0));

        CheckNormalization(GetNModes());  // E fields integrate to 1.0 for both TE and TM modes.

        if( aParam.has( "mode-maps" ) )
        {
        	PrintModeMaps(GetNModes(),1);
        }


/*
        if( aParam.has( "modemap-filename" ) )
        {

        }
*/


    	return true;
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
    	    		    	aField = fFieldCore->TE_E(GetDimR(), GetDimL(), l, m, n, r, theta, zKass,1);
    	    			}
    	    			else
    	    			{
    	    				aField = fFieldCore->TE_H(GetDimR(), GetDimL(), l, m, n, r, theta, zKass,1);
    	    			}
    	    		}
    	    		else
    	    		{
    	    			if (eField)
    	    			{
    	    				aField = fFieldCore->TM_E(GetDimR(), GetDimL(), l, m, n, r, theta, zKass,1);
    	    			}
    	    			else
    	    			{
    	    				aField = fFieldCore->TM_H(GetDimR(), GetDimL(), l, m, n, r, theta, zKass,1);
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

    std::vector<double> PozarCylindrical::TE_E(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool avgOverTheta)
    {

    	double z = zKass + L/2.;

    	// from Pozar
    	std::vector<double> TE_E;
    	double x_lm = GetBesselNKPrimeZeros(l,m);

    	double k1 = x_lm / R;
    	double k3 = n * LMCConst::Pi() / L;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );  // Pozar p. 291.
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * ( boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r) );
    	double tEr = 0.;
    	double tEtheta = 0.;

    	if ((!avgOverTheta)||(l==0))
    	{
        	tEr = -l * k/k1 * eta * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);
    		tEtheta = -k/k1 * eta * jPrime * cos(l*theta) * sin(k3*z);
    	}
    	else
    	{
        	tEr = -l * k/k1 * eta * jl_of_k1r_by_k1r * (2./LMCConst::Pi()) * sin(k3*z);
    		tEtheta = -k/k1 * eta * jPrime * (2./LMCConst::Pi()) * sin(k3*z);
    	}


    	TE_E.push_back(tEr);
    	TE_E.push_back(tEtheta);

        return TE_E;
    }

    std::vector<double> PozarCylindrical::TE_H(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool avgOverTheta)
    {

    	double z = zKass + L/2.;

    	// from Pozar
    	std::vector<double> TE_H;
    	double x_lm = GetBesselNKPrimeZeros(l,m);
    	double k1 = x_lm / R;
    	double k3 = n * LMCConst::Pi() / L;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * ( boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r) );
    	double tHz = boost::math::cyl_bessel_j(l, k1*r) * cos(l*theta) * sin(k3*z);
    	double tHr = 0.;
    	double tHtheta = 0.;

    	if ((!avgOverTheta)||(l==0))
    	{
        	tHr = -k3/k1 * jPrime * cos(l*theta) * cos(k3*z);
    		tHtheta = -l*k3/k1 * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z);
    	}
    	else
    	{
        	tHr = -k3/k1 * jPrime * (2./LMCConst::Pi()) * cos(k3*z);
    		tHtheta = -l*k3/k1 * jl_of_k1r_by_k1r * (2./LMCConst::Pi()) * cos(k3*z);
    	}

    	TE_H.push_back(tHr);  // r
    	TE_H.push_back(tHz);  // z
    	TE_H.push_back(tHtheta); // theta
    	return TE_H; // r, z, theta
    }

    std::vector<double> PozarCylindrical::TM_E(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool avgOverTheta)
    {
    	double z = zKass + L/2.;

    	// from Pozar
    	std::vector<double> TM_E;
    	double x_lm = GetBesselNKZeros(l,m);
    	double k1 = x_lm / R;
    	double k3 = n * LMCConst::Pi() / L;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );  // Pozar p. 291.
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * ( boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r) );
    	double tEz = eta * boost::math::cyl_bessel_j(l, k1*r) * cos(l*theta) * sin(k3*z);
    	double tEr = 0.;
    	double tEtheta = 0.;

    	if ((!avgOverTheta)||(l==0))
    	{
        	tEr = -k3/k1 * eta * jPrime * cos(l*theta) * cos(k3*z);
    		tEtheta = -l*k3/k1 * eta * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z);
    	}
    	else
    	{
        	tEr = -k3/k1 * eta * jPrime * (2./LMCConst::Pi()) * cos(k3*z);
    		tEtheta = -l*k3/k1 * eta * jl_of_k1r_by_k1r * (2./LMCConst::Pi()) * cos(k3*z);
    	}

    	TM_E.push_back(tEr); // r
    	TM_E.push_back(tEz);  // z
    	TM_E.push_back(tEtheta);  // theta
    	return TM_E; // r, z, theta
    }

    std::vector<double> PozarCylindrical::TM_H(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool avgOverTheta)
    {
    	double z = zKass + L/2.;

    	// from Pozar
    	std::vector<double> TM_H;
    	double x_lm = GetBesselNKZeros(l,m);
    	double k1 = x_lm / R;
    	double k3 = n * LMCConst::Pi() / L;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * ( boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r) );
    	double tHr = 0.;
    	double tHtheta = 0.;

    	if ((!avgOverTheta)||(l==0))
    	{
        	tHr = -l * k/k1  * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);
    		tHtheta = -k/k1 * jPrime * cos(l*theta) * sin(k3*z);
    	}
    	else
    	{
        	tHr = -l * k/k1  * jl_of_k1r_by_k1r * (2./LMCConst::Pi()) * sin(k3*z);
    		tHtheta = -k/k1 * jPrime * (2./LMCConst::Pi()) * sin(k3*z);
    	}

    	TM_H.push_back(tHr);  // r
    	TM_H.push_back(tHtheta); // theta
        return TM_H;
    }

    std::vector<double> CylindricalCavity::GetTE_E(int l, int m, int n, double r, double theta, double z, bool avgOverTheta)
    {
    	return fFieldCore->TE_E(GetDimR(),GetDimL(),l,m,n,r,0.,z,1);
    }

    std::vector<double> CylindricalCavity::GetNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP)
       {
       	double tR = tKassParticleXP[0];
       	double tZ = tKassParticleXP[2];
       	std::vector<double> tField;

       	tField = fFieldCore->TE_E(GetDimR(),GetDimL(),l,m,n,tR,0.,tZ,1);
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
    					getchar();
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



    void CylindricalCavity::PrintModeMaps(int nModes, bool bTE)
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
    					double r = (double)i/GetNPixels()*GetDimR();
    					for (unsigned j=0; j<GetNPixels()+1; j++)
    					{
    						double theta = (double)j/GetNPixels()*2.*LMCConst::Pi();
        					for (unsigned k=0; k<GetNPixels()+1; k++)
        					{
            				    double z = (double)k/GetNPixels()*GetDimL() - GetDimL()/2.;
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
    							fprintf(fp_E, "%10.4g %10.4g %10.4g %10.4g %10.4g\n", r, theta, z, tE.front()*normFactor, tE.back()*normFactor);
    							fprintf(fp_H, "%10.4g %10.4g %10.4g %10.4g %10.4g\n", r, theta, z, tH.front()*normFactor, tH.back()*normFactor);
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


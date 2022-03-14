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
    fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }
    CylindricalCavity::~CylindricalCavity() {}


    double CylindricalCavity::Integrate(int l, int m, int n, bool teMode, bool eField)
    {

    	std::vector<double> aField;
    	double r, theta, zPozar, zKass = 0.;
    	double dR = fInterface->fR/GetNPixels();
    	double dZ = fInterface->fL/GetNPixels();
    	double dTheta = 2.*LMCConst::Pi()/GetNPixels();
    	double tVolume = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<GetNPixels(); i++)
    		for (unsigned j=0; j<GetNPixels(); j++)
    			for (unsigned k=0; k<GetNPixels(); k++)
    			{
    	    		r = (double)i*dR;
    	    		theta = (double)j*dTheta;
    	    		zPozar = (double)k*dZ;
    	    		zKass = zPozar - fInterface->fL/2.;

    	    		if (teMode)
    	    		{
    	    			if (eField)
    	    			{
    	    		    	aField = TE_E(l, m, n, r, theta, zKass, GetCentralFrequency(),1);
    	    			}
    	    			else
    	    			{
    	    				aField = TE_H(l, m, n, r, theta, zKass, GetCentralFrequency(),1);
    	    			}
    	    		}
    	    		else
    	    		{
    	    			if (eField)
    	    			{
    	    				aField = TM_E(l, m, n, r, theta, zKass, GetCentralFrequency(),1);
    	    			}
    	    			else
    	    			{
    	    				aField = TM_H(l, m, n, r, theta, zKass, GetCentralFrequency(),1);
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

    double CylindricalCavity::GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP)
    {
    	// fix me (placeholder)
    	return 0.;
//    	return tKassParticleXP[7];  // fcyc
    }


    double CylindricalCavity::Z_TE(int l, int m, int n, double fcyc) const
    {
    	double Z_TE = 1.0;
    	double x_lm = fInterface->fBesselNKPrimeZeros[l][m];
    	double k1 = x_lm / fInterface->fR;
    	double k3 = n * LMCConst::Pi() / fInterface->fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double k0 = fcyc / LMCConst::C();
    	double A = 1.;
    	double B = k*k;
    	double C = k0*k0;
    	double Q = 1.;  // Q for suppressed modes.
    	if ((l==0)&&(m==1)&&(n==1)) Q = 1000.; // Q for lmn = 011.
    	double denom = (Q*B - Q*C - C) * (Q*B - Q*C - C) + C*C;
    	double real = (Q*A * (Q*B - Q*C - C)) / denom;
    	double imag = -Q*A*C / denom;

    	if ( k*k-k0*k0 != 0. )
    	{
    		// after Collin Foundations of M.E. Eq. 7.132
//    		Z_TE *= fcyc * LMCConst::MuNull() / (k0*k0 - k*k);  // neglect Q term.
    		Z_TE *= fcyc * LMCConst::MuNull() * sqrt( real*real + imag*imag ); // mag
    	}
    	return Z_TE;
    }

    double CylindricalCavity::Z_TM(int l, int m, int n, double fcyc) const
    {
    	double Z_TM = 1.0;
    	double x_lm = fInterface->fBesselNKZeros[l][m];
    	double k1 = x_lm / fInterface->fR;
    	double k3 = n * LMCConst::Pi() / fInterface->fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double k0 = fcyc / LMCConst::C();
    	double A = 1.;
    	double B = k*k;
    	double C = k0*k0;
    	double Q = 1.;  // Q for suppressed modes.
    	if ((l==0)&&(m==1)&&(n==1)) Q = 1000.; // Q for lmn = 011.
    	double denom = (Q*B - Q*C - C) * (Q*B - Q*C - C) + C*C;
    	double real = (Q*A * (Q*B - Q*C - C)) / denom;
    	double imag = -Q*A*C / denom;

    	if ( k*k-k0*k0 != 0. )
    	{
    		// after Collin Foundations of M.E. Eq. 7.132
//    		Z_TM *= fcyc * LMCConst::MuNull() / (k0*k0 - k*k);  // neglect Q term.
    		Z_TM *= fcyc * LMCConst::MuNull() * sqrt( real*real + imag*imag ); // mag
    	}

    	return Z_TM;
    }

    std::vector<double> CylindricalCavity::TE_E(int l, int m, int n, double r, double theta, double zKass, double fcyc, bool avgOverTheta) const
    {

    	double z = zKass + fInterface->fL/2.;

    	// from Pozar
    	std::vector<double> TE_E;
    	double x_lm = fInterface->fBesselNKPrimeZeros[l][m];
    	double k1 = x_lm / fInterface->fR;
    	double k3 = n * LMCConst::Pi() / fInterface->fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );  // Pozar p. 291.
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
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

    std::vector<double> CylindricalCavity::TE_H(int l, int m, int n, double r, double theta, double zKass, double fcyc, bool avgOverTheta) const
    {

    	double z = zKass + fInterface->fL/2.;

    	// from Pozar
    	std::vector<double> TE_H;
    	double x_lm = fInterface->fBesselNKPrimeZeros[l][m];
    	double k1 = x_lm / fInterface->fR;
    	double k3 = n * LMCConst::Pi() / fInterface->fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
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

    std::vector<double> CylindricalCavity::TM_E(int l, int m, int n, double r, double theta, double zKass, double fcyc, bool avgOverTheta) const
    {
    	double z = zKass + fInterface->fL/2.;

    	// from Pozar
    	std::vector<double> TM_E;
    	double x_lm = fInterface->fBesselNKZeros[l][m];
    	double k1 = x_lm / fInterface->fR;
    	double k3 = n * LMCConst::Pi() / fInterface->fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );  // Pozar p. 291.
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
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

    std::vector<double> CylindricalCavity::TM_H(int l, int m, int n, double r, double theta, double zKass, double fcyc, bool avgOverTheta) const
    {
    	double z = zKass + fInterface->fL/2.;

    	// from Pozar
    	std::vector<double> TM_H;
    	double x_lm = fInterface->fBesselNKZeros[l][m];
    	double k1 = x_lm / fInterface->fR;
    	double k3 = n * LMCConst::Pi() / fInterface->fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
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




} /* namespace locust */


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
    	double dR = fInterface->fR/fInterface->fnPixels;
    	double dZ = fInterface->fL/fInterface->fnPixels;
    	double dTheta = 2.*LMCConst::Pi()/fInterface->fnPixels;
    	double tVolume = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<fInterface->fnPixels; i++)
    		for (unsigned j=0; j<fInterface->fnPixels; j++)
    			for (unsigned k=0; k<fInterface->fnPixels; k++)
    			{
    	    		r = (double)i*dR;
    	    		theta = (double)j*dTheta;
    	    		zPozar = (double)k*dZ;
    	    		zKass = zPozar - fInterface->fL/2.;

    	    		if (teMode)
    	    		{
    	    			if (eField)
    	    			{
    	    		    	aField = TE_E(l, m, n, r, theta, zKass, GetCentralFrequency());
    	    			}
    	    			else
    	    			{
    	    				aField = TE_H(l, m, n, r, theta, zKass, GetCentralFrequency());
    	    			}
    	    		}
    	    		else
    	    		{
    	    			if (eField)
    	    			{
    	    				aField = TM_E(l, m, n, r, theta, zKass, GetCentralFrequency());
    	    			}
    	    			else
    	    			{
    	    				aField = TM_H(l, m, n, r, theta, zKass, GetCentralFrequency());
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

    std::vector<double> CylindricalCavity::TE_E(int l, int m, int n, double r, double theta, double zKass, double fcyc) const
    {

    	double z = zKass + fInterface->fL/2.;

    	// from Pozar
    	std::vector<double> TE_E;
    	double x_lm = fInterface->fBesselNKPrimeZeros[l][m];
    	double k1 = x_lm / fInterface->fR;
    	double k3 = n * LMCConst::Pi() / fInterface->fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double omega = LMCConst::C()*k;
    	double k0 = omega/LMCConst::C()*sqrt(LMCConst::MuNull()*LMCConst::EpsNull());
    	double eta = LMCConst::MuNull()*omega/LMCConst::C()/k0;  // Jackson 8.32
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double tEr = -l * k/k1 * eta * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
    	double tEtheta = -k/k1 * eta * jPrime * cos(l*theta) * sin(k3*z);
    	TE_E.push_back(tEr);
    	TE_E.push_back(tEtheta);
        return TE_E;
    }

    std::vector<double> CylindricalCavity::TE_H(int l, int m, int n, double r, double theta, double zKass, double fcyc) const
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
    	TE_H.push_back(-k3/k1 * jPrime * cos(l*theta) * cos(k3*z));
    	TE_H.push_back(-l*k3/k1 * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z));
    	TE_H.push_back(boost::math::cyl_bessel_j(l, k1*r) * cos(l*theta) * sin(k3*z));
        return TE_H;
    }

    std::vector<double> CylindricalCavity::TM_E(int l, int m, int n, double r, double theta, double zKass, double fcyc) const
    {
    	double z = zKass + fInterface->fL/2.;

    	// from Pozar
    	std::vector<double> TM_E;
    	double x_lm = fInterface->fBesselNKZeros[l][m];
    	double k1 = x_lm / fInterface->fR;
    	double k3 = n * LMCConst::Pi() / fInterface->fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double omega = LMCConst::C()*k;
    	double k0 = omega/LMCConst::C()*sqrt(LMCConst::MuNull()*LMCConst::EpsNull());
    	double eta = LMCConst::C()/LMCConst::EpsNull()/omega*k0;  // Jackson 8.32
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
    	TM_E.push_back(-k3/k1 * eta * jPrime * cos(l*theta) * cos(k3*z));
    	TM_E.push_back(-l*k3/k1 * eta * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z));
    	TM_E.push_back(eta * boost::math::cyl_bessel_j(l, k1*r) * cos(l*theta) * sin(k3*z));
        return TM_E;
    }

    std::vector<double> CylindricalCavity::TM_H(int l, int m, int n, double r, double theta, double zKass, double fcyc) const
    {
    	double z = zKass + fInterface->fL/2.;

    	// from Pozar
    	std::vector<double> TM_H;
    	double x_lm = fInterface->fBesselNKZeros[l][m];
    	double k1 = x_lm / fInterface->fR;
    	double k3 = n * LMCConst::Pi() / fInterface->fL;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double tHr = -l * k/k1  * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);
    	double jPrime = 1./2. * boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r);
    	double tHtheta = -k/k1 * jPrime * cos(l*theta) * sin(k3*z);
    	TM_H.push_back(tHr);
    	TM_H.push_back(tHtheta);
        return TM_H;
    }




} /* namespace locust */


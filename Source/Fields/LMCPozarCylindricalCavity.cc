/*
 * LMCPozarCylindricalCavity.cc
 *
 *  Created on: Jun 9, 2021
 *      Author: pslocum
 */

#include "LMCPozarCylindricalCavity.hh"


namespace locust
{

    LOGGER( lmclog, "PozarCylindricalCavity" );
    PozarCylindricalCavity::PozarCylindricalCavity()
    {
    }

    PozarCylindricalCavity::~PozarCylindricalCavity(){}



    std::vector<double> PozarCylindricalCavity::TE_E(double R, double twoPi, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {

    	double z = zKass + L/2.;

    	if ((r > R) || (fabs(zKass) > L/2.))  // outside the cavity
    	{
    		return {0.,0.};
    	}
    	else
    	{

    	// from Pozar
    	std::vector<double> TE_E;
    	double x_lm = GetBesselNKPrimeZeros(l,m);

    	double k1 = x_lm / R;
    	double k3 = n * LMCConst::Pi() / L;
    	double k = 1.0; //pow(k1*k1+k3*k3,0.5);
    	double eta = 1.0; //sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );  // Pozar p. 291.
    	double jl_of_k1r_by_k1r = 1./2. * (boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r));

    	double jPrime = jl_of_k1r_by_k1r;
    	double tEr = l * k/k1/k1 * eta * boost::math::cyl_bessel_j(l, k1*r) * sin(l*theta) * sin(k3*z) / (r + 1e-17);
    	double tEtheta = k/k1 * eta * jPrime * cos(l*theta) * sin(k3*z);

    	if ((includeOtherPols)&&(l>0))
    	{
    		//modifies both r and theta components of TE field.
    		double dTheta = LMCConst::Pi() / 2.0 / (double)l;
    		std::vector<double> tPolarization = this->TE_E(R,2.*LMCConst::Pi(),L,l,m,n,r,theta+dTheta,zKass,0);
    		tEr = tEr*sin((double)l*theta) + tPolarization[0]*cos((double)l*theta) ;
    		tEtheta = tEtheta*sin((double)l*(theta+dTheta)) + tPolarization[1]*cos((double)l*(theta+dTheta)) ;
    	}

		double norm = R * sqrt(LMCConst::Pi() * L) / 2 / k1;
		norm *= sqrt(1 - l*l / x_lm / x_lm) * boost::math::cyl_bessel_j(l, x_lm);

		if (l == 0)
		{
			norm *= sqrt(2);
		}

    	TE_E.push_back(tEr / norm);
    	TE_E.push_back(tEtheta / norm);
		
        return TE_E;
    	}
    }

    std::vector<double> PozarCylindricalCavity::TE_H(double R, double twoPi, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {

    	double z = zKass + L/2.;

    	if ((r > R) || (fabs(zKass) > L/2.))  // outside the cavity
    	{
    		return {0.,0.,0.};
    	}
    	else
    	{

    	// from Pozar
    	std::vector<double> TE_H;
    	double x_lm = GetBesselNKPrimeZeros(l,m);
    	double k1 = x_lm / R;
    	double k3 = n * LMCConst::Pi() / L;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * ( boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r) );
    	double tHz = boost::math::cyl_bessel_j(l, k1*r) * cos(l*theta) * sin(k3*z);
    	double tHr = -k3/k1 * jPrime * cos(l*theta) * cos(k3*z);
    	double tHtheta = -l*k3/k1 * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z);

    	if ((includeOtherPols)&&(l>0))
    	{
    		double dTheta = LMCConst::Pi() / 2.0 / (double)l;
    		std::vector<double> tPolarization = this->TE_H(R,2.*LMCConst::Pi(),L,l,m,n,r,theta+dTheta,zKass,0);
    		tHr = tHr*sin((double)l*(theta+dTheta)) + tPolarization[0]*cos((double)l*(theta+dTheta)) ;
    		tHz = tHz*sin((double)l*(theta+dTheta)) + tPolarization[1]*cos((double)l*(theta+dTheta)) ;
    		tHtheta = tHtheta*sin((double)l*theta) + tPolarization[2]*cos((double)l*theta) ;
    	}

    	TE_H.push_back(tHr);  // r
    	TE_H.push_back(tHz);  // z
    	TE_H.push_back(tHtheta); // theta
    	return TE_H; // r, z, theta
    	}
    }

    std::vector<double> PozarCylindricalCavity::TM_E(double R, double twoPi, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
    	double z = zKass + L/2.;

    	if ((r > R) || (fabs(zKass) > L/2.))  // outside the cavity
    	{
    		return {0.,0.,0.};
    	}
    	else
    	{

    	// from Pozar
    	std::vector<double> TM_E;
    	double x_lm = GetBesselNKZeros(l,m);
    	double k1 = x_lm / R;
    	double k3 = n * LMCConst::Pi() / L;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );  // Pozar p. 291.
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * ( boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r) );
    	double tEz = eta * boost::math::cyl_bessel_j(l, k1*r) * cos(l*theta) * cos(k3*z);
    	double tEr = -k3/k1 * eta * jPrime * cos(l*theta) * sin(k3*z);
    	double tEtheta = -l*k3/k1 * eta * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);

    	if ((includeOtherPols)&&(l>0))
    	{
                double dTheta = LMCConst::Pi() / 2.0 / (double)l;
                std::vector<double> tPolarization = this->TM_E(R,2.*LMCConst::Pi(),L,l,m,n,r,theta+dTheta,zKass,0);
                tEr = tEr*sin((double)l*(theta+dTheta)) + tPolarization[0]*cos((double)l*(theta+dTheta)) ;
                tEz = tEz*sin((double)l*(theta+dTheta)) + tPolarization[1]*cos((double)l*(theta+dTheta)) ;
                tEtheta = tEtheta*sin((double)l*theta) + tPolarization[2]*cos((double)l*theta) ;
    	}

    	TM_E.push_back(tEr); // r
    	TM_E.push_back(tEz);  // z
    	TM_E.push_back(tEtheta);  // theta
    	return TM_E; // r, z, theta
    	}
    }

    std::vector<double> PozarCylindricalCavity::TM_H(double R, double twoPi, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
    {
    	double z = zKass + L/2.;

    	if ((r > R) || (fabs(zKass) > L/2.))  // outside the cavity
    	{
    		return {0.,0.};
    	}
    	else
    	{

    	// from Pozar
    	std::vector<double> TM_H;
    	double x_lm = GetBesselNKZeros(l,m);
    	double k1 = x_lm / R;
    	double k3 = n * LMCConst::Pi() / L;
    	double k = pow(k1*k1+k3*k3,0.5);
    	double jl_of_k1r_by_k1r = 1./(2.*l) * (boost::math::cyl_bessel_j(l-1, k1*r) + boost::math::cyl_bessel_j(l+1, k1*r));
    	double jPrime = 1./2. * ( boost::math::cyl_bessel_j(l-1, k1*r) - boost::math::cyl_bessel_j(l+1, k1*r) );
    	double tHr = -l * k/k1  * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z);
    	double tHtheta = -k/k1 * jPrime * cos(l*theta) * cos(k3*z);

    	if ((includeOtherPols)&&(l>0))
    	{
                //modifies both r and theta components of TM field.
                double dTheta = LMCConst::Pi() / 2.0 / (double)l;
                std::vector<double> tPolarization = this->TM_H(R,2.*LMCConst::Pi(),L,l,m,n,r,theta+dTheta,zKass,0);
                tHr = tHr*sin((double)l*theta) + tPolarization[0]*cos((double)l*theta) ;
                tHtheta = tHtheta*sin((double)l*(theta+dTheta)) + tPolarization[1]*cos((double)l*(theta+dTheta)) ;
    	}

    	TM_H.push_back(tHr);  // r
    	TM_H.push_back(tHtheta); // theta
        return TM_H;
    	}
    }

} /* namespace locust */


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



    std::vector<double> PozarCylindricalCavity::TE_E(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
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

    	tEr = -l * k/k1 * eta * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);
    	tEtheta = -k/k1 * eta * jPrime * cos(l*theta) * sin(k3*z);

    	if ((includeOtherPols)&&(l>0))
    	{
    		//modifies both r and phi components of TE field.
    		double dPhi = LMCConst::Pi() / 2.0 / (double)l;
    		std::vector<double> tPolarization{tEr, tEtheta};
    		tEr = tEr*sin((double)l*theta) + tPolarization[0]*cos((double)l*theta) ;
    		tEtheta = tEtheta*sin((double)l*(theta+dPhi)) + tPolarization[1]*cos((double)l*(theta+dPhi)) ;
    	}

    	TE_E.push_back(tEr);
    	TE_E.push_back(tEtheta);

        return TE_E;
    }

    std::vector<double> PozarCylindricalCavity::TE_H(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
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

    	if ((!includeOtherPols)||(l==0))
    	{
        	tHr = -k3/k1 * jPrime * cos(l*theta) * cos(k3*z);
    		tHtheta = -l*k3/k1 * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z);
    	}
    	else
    	{
    		LERROR(lmclog,"This superposition has not yet been implemented.");
    		exit(-1);
    		// Possible suggestion:
    		// Here we can implement the superposition with other polarities of the same mode.
    		// The superposition can be done either in this function itself, or with some kind
    		// of new helper function in this class.
    	}

    	TE_H.push_back(tHr);  // r
    	TE_H.push_back(tHz);  // z
    	TE_H.push_back(tHtheta); // theta
    	return TE_H; // r, z, theta
    }

    std::vector<double> PozarCylindricalCavity::TM_E(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
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

    	if ((!includeOtherPols)||(l==0))
    	{
        	tEr = -k3/k1 * eta * jPrime * cos(l*theta) * cos(k3*z);
    		tEtheta = -l*k3/k1 * eta * jl_of_k1r_by_k1r * sin(l*theta) * cos(k3*z);
    	}
    	else
    	{
    		LERROR(lmclog,"This superposition has not yet been implemented.");
    		exit(-1);
    		// Possible suggestion:
    		// Here we can implement the superposition with other polarities of the same mode.
    		// The superposition can be done either in this function itself, or with some kind
    		// of new helper function in this class.
    	}

    	TM_E.push_back(tEr); // r
    	TM_E.push_back(tEz);  // z
    	TM_E.push_back(tEtheta);  // theta
    	return TM_E; // r, z, theta
    }

    std::vector<double> PozarCylindricalCavity::TM_H(double R, double L, int l, int m, int n, double r, double theta, double zKass, bool includeOtherPols)
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

    	if ((!includeOtherPols)||(l==0))
    	{
        	tHr = -l * k/k1  * jl_of_k1r_by_k1r * sin(l*theta) * sin(k3*z);
    		tHtheta = -k/k1 * jPrime * cos(l*theta) * sin(k3*z);
    	}
    	else
    	{
    		LERROR(lmclog,"This superposition has not yet been implemented.");
    		exit(-1);
    		// Possible suggestion:
    		// Here we can implement the superposition with other polarities of the same mode.
    		// The superposition can be done either in this function itself, or with some kind
    		// of new helper function in this class.
    	}

    	TM_H.push_back(tHr);  // r
    	TM_H.push_back(tHtheta); // theta
        return TM_H;
    }

} /* namespace locust */


/*
 * LMCPozarRectangularCavity.cc
 *
 *  Created on: Nov 16, 2023
 *      Author: pslocum
 */

#include "LMCPozarRectangularCavity.hh"


namespace locust
{
    LOGGER( lmclog, "PozarRectangularCavity" );
    PozarRectangularCavity::PozarRectangularCavity()
    {
    }
    PozarRectangularCavity::~PozarRectangularCavity() {}


    std::vector<double> PozarRectangularCavity::TE_E(double dimX, double dimY, double dimZ, int l, int m, int n, double xKass, double yKass, double zKass, bool includeOtherPols)
    {

    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;
    	double z = zKass + dimZ/2.;

    	if ((fabs(xKass) > dimX/2.) || (fabs(yKass) > dimY/2.) || (fabs(zKass) > dimZ/2.))
    	{
    		return {0.,0.};
    	}
    	else
    	{

    	// from Pozar
    	std::vector<double> TE_E;

    	double k1 = l * LMCConst::Pi() / dimX;
    	double k2 = m * LMCConst::Pi() / dimY;
    	double k3 = n * LMCConst::Pi() / dimZ;
    	double kc = pow(k1*k1+k2*k2+k3*k3,0.5);
    	double beta = k3;
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );  // Pozar p. 113
    	double Z_TE = kc * eta / beta;

    	double tEx = Z_TE * beta * k2 / kc / kc * cos(k1*x) * sin(k2*y) * sin(k3*z);
    	double tEy = -Z_TE * beta * k1 / kc / kc * sin(k1*x) * cos(k2*y) * sin(k3*z);
    	double tEz = 0.; // TE mode.

    	TE_E.push_back(tEx);
    	TE_E.push_back(tEy);

        return TE_E;
    	}
    }


    std::vector<double> PozarRectangularCavity::TE_H(double dimX, double dimY, double dimZ, int l, int m, int n, double xKass, double yKass, double zKass, bool includeOtherPols)
    {
    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;
    	double z = zKass + dimZ/2.;

    	if ((fabs(xKass) > dimX/2.) || (fabs(yKass) > dimY/2.) || (fabs(zKass) > dimZ/2.))
    	{
    		return {0.,0.};
    	}
    	else
    	{

    	// from Pozar
    	std::vector<double> TE_H;

    	double k1 = l * LMCConst::Pi() / dimX;
    	double k2 = m * LMCConst::Pi() / dimY;
    	double k3 = n * LMCConst::Pi() / dimZ;
    	double kc = pow(k1*k1+k2*k2+k3*k3,0.5);
    	double beta = k3;
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );  // Pozar p. 113
    	double Z_TE = kc * eta / beta;

    	double tHx = beta * k1 / kc / kc * sin(k1*x) * cos(k2*y) * cos(k3*z);
    	double tHy = beta * k2 / kc / kc * cos(k1*x) * sin(k2*y) * cos(k3*z);
    	double tHz = beta * (k1*k1 + k2*k2) / kc / kc * Z_TE / kc / eta  * cos(k1*x) * cos(k2*y) * sin(k3*z); // 6.42c

    	TE_H.push_back(tHx);
    	TE_H.push_back(tHz);
    	TE_H.push_back(tHy);

        return TE_H;
    	}
    }


    std::vector<double> PozarRectangularCavity::TM_E(double dimX, double dimY, double dimZ, int l, int m, int n, double xKass, double yKass, double zKass, bool includeOtherPols)
    {
    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;
    	double z = zKass + dimZ/2.;

    	if ((fabs(xKass) > dimX/2.) || (fabs(yKass) > dimY/2.) || (fabs(zKass) > dimZ/2.))
    	{
    		return {0.,0.};
    	}
    	else
    	{

    	// from Pozar
    	std::vector<double> TM_E;

       	double k1 = l * LMCConst::Pi() / dimX;
        double k2 = m * LMCConst::Pi() / dimY;
        double k3 = n * LMCConst::Pi() / dimZ;
        double kc = pow(k1*k1+k2*k2+k3*k3,0.5);
    	double beta = k3;
        double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );  // Pozar p. 113
    	double Z_TE = kc * eta / beta;

    	double tEx = -beta * k1 / kc / kc * cos(k1*x) * sin(k2*y) * sin(k3*z);
    	double tEy = -beta * k2 / kc / kc * sin(k1*x) * cos(k2*y) * sin(k3*z);
    	double tEz =  Z_TE * beta * (k1*k1 + k2*k2) / kc / kc / kc / eta * sin(k1*x) * sin(k2*y) * cos(k3*z);

        TM_E.push_back(tEx);
        TM_E.push_back(tEz);
        TM_E.push_back(tEy);

    	return TM_E;
    	}
    }

    std::vector<double> PozarRectangularCavity::TM_H(double dimX, double dimY, double dimZ, int l, int m, int n, double xKass, double yKass, double zKass, bool includeOtherPols)
    {
    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;
    	double z = zKass + dimZ/2.;

    	if ((fabs(xKass) > dimX/2.) || (fabs(yKass) > dimY/2.) || (fabs(zKass) > dimZ/2.))
    	{
    		return {0.,0.};
    	}
    	else
    	{

    	// from Pozar
    	std::vector<double> TM_H;

    	double k1 = l * LMCConst::Pi() / dimX;
    	double k2 = m * LMCConst::Pi() / dimY;
    	double k3 = n * LMCConst::Pi() / dimZ;
    	double kc = pow(k1*k1+k2*k2+k3*k3,0.5);
    	double beta = k3;
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );
    	double Z_TM = beta * eta / kc;

    	double tHx = beta * k2 / kc / kc / Z_TM * sin(k1*x) * cos(k2*y) * cos(k3*z);
    	double tHy = -beta * k1 / kc / kc / Z_TM * cos(k1*x) * sin(k2*y) * cos(k3*z);
    	double tHz = 0.; // TM mode.

    	TM_H.push_back(tHx);
    	TM_H.push_back(tHy);

        return TM_H;
    	}
    }


} /* namespace locust */


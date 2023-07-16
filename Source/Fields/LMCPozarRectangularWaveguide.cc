/*
 * LMCPozarRectangularWaveguide.cc
 *
 *  Created on: Apr 11, 2023
 *      Author: pslocum
 */

#include "LMCPozarRectangularWaveguide.hh"


namespace locust
{
    LOGGER( lmclog, "PozarRectangularWaveguide" );
    PozarRectangularWaveguide::PozarRectangularWaveguide()
    {
    }
    PozarRectangularWaveguide::~PozarRectangularWaveguide() {}


    std::vector<double> PozarRectangularWaveguide::TE_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc)
    {

    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;

    	if ((fabs(x) > dimX/2.) || (fabs(y) > dimY/2.))
    	{
    		return {0.,0.};
    	}
    	else
    	{

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
    }


    std::vector<double> PozarRectangularWaveguide::TE_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc)
    {
    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;

    	if ((fabs(x) > dimX/2.) || (fabs(y) > dimY/2.))
    	{
    		return {0.,0.};
    	}
    	else
    	{

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
    }


    std::vector<double> PozarRectangularWaveguide::TM_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc)
    {
    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;

    	if ((fabs(x) > dimX/2.) || (fabs(y) > dimY/2.))
    	{
    		return {0.,0.};
    	}
    	else
    	{

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
    }

    std::vector<double> PozarRectangularWaveguide::TM_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc)
    {
    	double x = xKass + dimX/2.;
    	double y = yKass + dimY/2.;

    	if ((fabs(x) > dimX/2.) || (fabs(y) > dimY/2.))
    	{
    		return {0.,0.};
    	}
    	else
    	{

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
    }


} /* namespace locust */


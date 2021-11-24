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


    double RectangularWaveguide::Integrate(int l, int m, int n, bool teMode, bool eField)
    {

    	std::vector<double> aField;
    	double xPozar, yPozar, xKass, yKass = 0.;
    	double dX = fInterface->fX/fInterface->fnPixels;
    	double dY = fInterface->fY/fInterface->fnPixels;
    	double tArea = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<fInterface->fnPixels; i++)
    		for (unsigned j=0; j<fInterface->fnPixels; j++)
    		{
    	    	xPozar = (double)i*dX;
    	    	yPozar = (double)j*dY;
   	    		xKass = xPozar - fInterface->fX/2.;
   	    		yKass = yPozar - fInterface->fY/2.;

    	    	if (teMode)
    	    	{
    	    		if (eField)
    	    		{
    	    		    aField = TE_E(m, n, xKass, yKass, GetCentralFrequency());
    	    		}
    	    		else
    	    		{
    	    			aField = TE_H(m, n, xKass, yKass, GetCentralFrequency());
    	    		}
    	    	}
    	    	else
    	    	{
    	    		if (eField)
    	    		{
    	    			aField = TM_E(m, n, xKass, yKass, GetCentralFrequency());
    	    		}
    	    		else
    	    		{
    	    			aField = TM_H(m, n, xKass, yKass, GetCentralFrequency());
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

    			tIntegral += aFieldMagSq*dX*dY;
//    		    tArea += dX*dY;  // sanity check area integral.
    		}
//    	printf("tArea is %g\n", tArea); getchar();
    	return tIntegral;
    }

    std::vector<double> RectangularWaveguide::TE_E(int m, int n, double xKass, double yKass, double fcyc) const
    {

    	double x = xKass + fInterface->fX/2.;
    	double y = yKass + fInterface->fY/2.;

    	// from Pozar
    	std::vector<double> TE_E;
    	double k1 = m * LMCConst::Pi() / fInterface->fX;
    	double k2 = n * LMCConst::Pi() / fInterface->fY;
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double omega = LMCConst::C()*kc;  // TO-DO:  Make this frequency-dependent.
    	double k = omega / sqrt(LMCConst::EpsNull()*LMCConst::MuNull());

    	double tEx = omega*LMCConst::MuNull()*n*LMCConst::Pi()/kc/kc/fInterface->fY * cos(k1*x) * sin(k2*y);
    	double tEy = -omega*LMCConst::MuNull()*m*LMCConst::Pi()/kc/kc/fInterface->fX * sin(k1*x) * cos(k2*y);
    	TE_E.push_back(tEx);
    	TE_E.push_back(tEy);
        return TE_E;
    }

    std::vector<double> RectangularWaveguide::TE_H(int m, int n, double xKass, double yKass, double fcyc) const
    {
    	double x = xKass + fInterface->fX/2.;
    	double y = yKass + fInterface->fY/2.;

    	// from Pozar
    	std::vector<double> TE_H;
    	double k1 = m * LMCConst::Pi() / fInterface->fX;
    	double k2 = n * LMCConst::Pi() / fInterface->fY;
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double omega = LMCConst::C()*kc;  // TO-DO:  Make this frequency-dependent.
    	double k = omega / sqrt(LMCConst::EpsNull()*LMCConst::MuNull());
    	double beta = sqrt(k*k - kc*kc);

    	double tHx = beta*m*LMCConst::Pi()/kc/kc/fInterface->fX * sin(k1*x) * cos(k2*y);
    	double tHy = beta*n*LMCConst::Pi()/kc/kc/fInterface->fY * cos(k1*x) * sin(k2*y);
    	TE_H.push_back(tHx);
    	TE_H.push_back(tHy);
        return TE_H;
    }

    std::vector<double> RectangularWaveguide::TM_E(int m, int n, double xKass, double yKass, double fcyc) const
    {
    	double x = xKass + fInterface->fX/2.;
    	double y = yKass + fInterface->fY/2.;

    	// from Pozar
    	std::vector<double> TM_E;
    	double k1 = m * LMCConst::Pi() / fInterface->fX;
    	double k2 = n * LMCConst::Pi() / fInterface->fY;
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double omega = LMCConst::C()*kc;  // TO-DO:  Make this frequency-dependent.
    	double k = omega / sqrt(LMCConst::EpsNull()*LMCConst::MuNull());
    	double beta = sqrt(k*k - kc*kc);

    	double tEx = beta*m*LMCConst::Pi()/kc/kc/fInterface->fX * cos(k1*x) * sin(k2*y);
    	double tEy = beta*n*LMCConst::Pi()/kc/kc/fInterface->fY * sin(k1*x) * cos(k2*y);
    	TM_E.push_back(tEx);
    	TM_E.push_back(tEy);
        return TM_E;
    }

    std::vector<double> RectangularWaveguide::TM_H(int m, int n, double xKass, double yKass, double fcyc) const
    {
    	double x = xKass + fInterface->fX/2.;
    	double y = yKass + fInterface->fY/2.;

    	// from Pozar
    	std::vector<double> TM_H;
    	double k1 = m * LMCConst::Pi() / fInterface->fX;
    	double k2 = n * LMCConst::Pi() / fInterface->fY;
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double omega = LMCConst::C()*kc;  // TO-DO:  Make this frequency-dependent.
    	double k = omega / sqrt(LMCConst::EpsNull()*LMCConst::MuNull());

    	double tHx = omega*LMCConst::EpsNull()*n*LMCConst::Pi()/kc/kc/fInterface->fY * sin(k1*x) * cos(k2*y);
    	double tHy = omega*LMCConst::EpsNull()*m*LMCConst::Pi()/kc/kc/fInterface->fX * cos(k1*x) * sin(k2*y);
    	TM_H.push_back(tHx);
    	TM_H.push_back(tHy);
        return TM_H;
    }




} /* namespace locust */


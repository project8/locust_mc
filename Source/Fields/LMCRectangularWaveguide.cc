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
    	double dX = fInterface->fX/GetNPixels();
    	double dY = fInterface->fY/GetNPixels();
    	double tArea = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<GetNPixels(); i++)
    		for (unsigned j=0; j<GetNPixels(); j++)
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


    double RectangularWaveguide::GetGroupVelocity(int m, int n, double fcyc)
    {
    	double CutOffFrequency = 0.;
    	if ((m<2)&&(n<1))  // most likely case
    	{
    		// rad/s
    		CutOffFrequency = LMCConst::C() * LMCConst::Pi() / fInterface->fX;
    	}
    	else  // general case
    	{
    		// rad/s
    		CutOffFrequency = LMCConst::C() *
    				sqrt(pow(m*LMCConst::Pi()/fInterface->fX,2.) + sqrt(pow(n*LMCConst::Pi()/fInterface->fY,2.)));
    	}
        double GroupVelocity = LMCConst::C() * pow( 1. - pow(CutOffFrequency/fcyc, 2.) , 0.5);
        //        printf("GroupVelocity is %g\n", GroupVelocity); getchar();
        return GroupVelocity;
    }


    double RectangularWaveguide::GetDopplerFrequency(int l, int m, int n, std::vector<double> tKassParticleXP)
    {
    	double fcyc = tKassParticleXP[7];
    	double groupVelocity = GetGroupVelocity(m,n,fcyc);
    	double zVelocity = tKassParticleXP[5];
        double gammaZ = 1.0 / pow(1.0-pow(zVelocity/groupVelocity,2.),0.5);
        double fPrime = fcyc * gammaZ * (1.+zVelocity/groupVelocity);
    	return fPrime;
    }


    double RectangularWaveguide::Z_TE(int l, int m, int n, double fcyc) const
    {
    	double k1 = m * LMCConst::Pi() / fInterface->fX;
    	double k2 = n * LMCConst::Pi() / fInterface->fY;
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );
    	double k = fcyc / LMCConst::C();
    	double beta = sqrt(k*k - kc*kc);

    	double Z_TE = k*eta/beta;  // This is 448 ohms for TE10 at 25.9 GHz.
    	return 2. * LMCConst::Pi() * Z_TE / LMCConst::C() / 1.e2; // Jackson Eq. 8.140, 1.e2 is m/s -> cm/s
    }

    double RectangularWaveguide::Z_TM(int l, int m, int n, double fcyc) const
    {
    	double k1 = m * LMCConst::Pi() / fInterface->fX;
    	double k2 = n * LMCConst::Pi() / fInterface->fY;
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double eta = sqrt( LMCConst::MuNull() / LMCConst::EpsNull() );
    	double k = fcyc / LMCConst::C();
    	double beta = sqrt(k*k - kc*kc);

    	double Z_TM = beta*eta/k;
    	return 2. * LMCConst::Pi() * Z_TM / LMCConst::C() / 1.e2; // Jackson Eq. 8.140, 1.e2 is m/s -> cm/s
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

    	double tEx = fcyc*LMCConst::MuNull()*n*LMCConst::Pi()/kc/kc/fInterface->fY * cos(k1*x) * sin(k2*y);
    	double tEy = -fcyc*LMCConst::MuNull()*m*LMCConst::Pi()/kc/kc/fInterface->fX * sin(k1*x) * cos(k2*y);
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
    	double k = fcyc * sqrt(LMCConst::EpsNull()*LMCConst::MuNull());
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
    	double k = fcyc * sqrt(LMCConst::EpsNull()*LMCConst::MuNull());
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

    	double tHx = fcyc*LMCConst::EpsNull()*n*LMCConst::Pi()/kc/kc/fInterface->fY * sin(k1*x) * cos(k2*y);
    	double tHy = fcyc*LMCConst::EpsNull()*m*LMCConst::Pi()/kc/kc/fInterface->fX * cos(k1*x) * sin(k2*y);
    	TM_H.push_back(tHx);
    	TM_H.push_back(tHy);
        return TM_H;
    }




} /* namespace locust */


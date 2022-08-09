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
    	double dX = GetDimX()/GetNPixels();
    	double dY = GetDimY()/GetNPixels();
    	double tArea = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<GetNPixels(); i++)
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


    std::vector<double> RectangularWaveguide::TE_E(int m, int n, double xKass, double yKass, double fcyc) const
    {

    	double x = xKass + GetDimX()/2.;
    	double y = yKass + GetDimY()/2.;

    	// from Pozar
    	std::vector<double> TE_E;
    	double k1 = m * LMCConst::Pi() / GetDimX();
    	double k2 = n * LMCConst::Pi() / GetDimY();
    	double kc = pow(k1*k1+k2*k2,0.5);

    	double tEx = fcyc*LMCConst::MuNull()*n*LMCConst::Pi()/kc/kc/GetDimY() * cos(k1*x) * sin(k2*y);
    	double tEy = -fcyc*LMCConst::MuNull()*m*LMCConst::Pi()/kc/kc/GetDimX() * sin(k1*x) * cos(k2*y);
    	TE_E.push_back(tEx);
    	TE_E.push_back(tEy);
        return TE_E;
    }

    std::vector<double> RectangularWaveguide::TE_H(int m, int n, double xKass, double yKass, double fcyc) const
    {
    	double x = xKass + GetDimX()/2.;
    	double y = yKass + GetDimY()/2.;

    	// from Pozar
    	std::vector<double> TE_H;
    	double k1 = m * LMCConst::Pi() / GetDimX();
    	double k2 = n * LMCConst::Pi() / GetDimY();
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double k = fcyc * sqrt(LMCConst::EpsNull()*LMCConst::MuNull());
    	double beta = sqrt(k*k - kc*kc);

    	double tHx = beta*m*LMCConst::Pi()/kc/kc/GetDimX() * sin(k1*x) * cos(k2*y);
    	double tHy = beta*n*LMCConst::Pi()/kc/kc/GetDimY() * cos(k1*x) * sin(k2*y);

    	TE_H.push_back(tHx);
    	TE_H.push_back(tHy);
        return TE_H;
    }


    std::vector<double> RectangularWaveguide::TM_E(int m, int n, double xKass, double yKass, double fcyc) const
    {
    	double x = xKass + GetDimX()/2.;
    	double y = yKass + GetDimY()/2.;

    	// from Pozar
    	std::vector<double> TM_E;
    	double k1 = m * LMCConst::Pi() / GetDimX();
    	double k2 = n * LMCConst::Pi() / GetDimY();
    	double kc = pow(k1*k1+k2*k2,0.5);
    	double k = fcyc * sqrt(LMCConst::EpsNull()*LMCConst::MuNull());
    	double beta = sqrt(k*k - kc*kc);

    	double tEx = beta*m*LMCConst::Pi()/kc/kc/GetDimX() * cos(k1*x) * sin(k2*y);
    	double tEy = beta*n*LMCConst::Pi()/kc/kc/GetDimY() * sin(k1*x) * cos(k2*y);
    	TM_E.push_back(tEx);
    	TM_E.push_back(tEy);
        return TM_E;
    }

    std::vector<double> RectangularWaveguide::TM_H(int m, int n, double xKass, double yKass, double fcyc) const
    {
    	double x = xKass + GetDimX()/2.;
    	double y = yKass + GetDimY()/2.;

    	// from Pozar
    	std::vector<double> TM_H;
    	double k1 = m * LMCConst::Pi() / GetDimX();
    	double k2 = n * LMCConst::Pi() / GetDimY();
    	double kc = pow(k1*k1+k2*k2,0.5);

    	double tHx = fcyc*LMCConst::EpsNull()*n*LMCConst::Pi()/kc/kc/GetDimY() * sin(k1*x) * cos(k2*y);
    	double tHy = fcyc*LMCConst::EpsNull()*m*LMCConst::Pi()/kc/kc/GetDimX() * cos(k1*x) * sin(k2*y);
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
      	std::vector<double> tTE_E_electron = this->TE_E(m,n,tX,tY,fcyc);
  		double normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];

  		auto it = tTE_E_electron.begin();
  		while (it != tTE_E_electron.end())
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




} /* namespace locust */


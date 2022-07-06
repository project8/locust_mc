
#include "LMCFieldCalculator.hh"
#include "KSParticleFactory.h"
#include <algorithm>

#include "KSInteractionsMessage.h"
#include <limits>

using std::numeric_limits;

using namespace Kassiopeia;
namespace locust
{

    FieldCalculator::FieldCalculator() :
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }
    FieldCalculator::FieldCalculator( const FieldCalculator& aCopy ) :
            fInterface( aCopy.fInterface )
    {
    }
    FieldCalculator* FieldCalculator::Clone() const
    {
        return new FieldCalculator( *this );
    }
    FieldCalculator::~FieldCalculator()
    {
    }


    double FieldCalculator::GetGroupVelocityTM01(Kassiopeia::KSParticle& aFinalParticle)
    {
        const double SpeedOfLight = LMCConst::C();
        double CutOffFrequency = 2. * LMCConst::Pi() * SpeedOfLight * 2.405 / 2. / LMCConst::Pi() / 0.00502920; // rad/s
        double cyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * sqrt( 1. - pow(CutOffFrequency/(2. * LMCConst::Pi() * cyclotronFrequency), 2.) );
        return GroupVelocity;
    }


    double FieldCalculator::GetGroupVelocityTE10(Kassiopeia::KSParticle& aFinalParticle)  // Phase 1
    {
        double SpeedOfLight = LMCConst::C(); // m/s
        double CutOffFrequency = SpeedOfLight * LMCConst::Pi() / fInterface->fField->GetDimX(); // a in m
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * pow( 1. - pow(CutOffFrequency/(2.*LMCConst::Pi()*fcyc), 2.) , 0.5);
        //        printf("GroupVelocity is %g\n", GroupVelocity); getchar();
        return GroupVelocity;
    }


    double FieldCalculator::GetCouplingFactorTM01(Kassiopeia::KSParticle& aFinalParticle)
    {
        double kc = 2.405/0.00502920;
        double x = aFinalParticle.GetPosition().GetX();
        double y = aFinalParticle.GetPosition().GetY();
        double r = sqrt(x*x + y*y);

    // power fraction. 48.1 is numerical normalization
   // of J \cdot E after time averaging as in Collin IEEE paper.
   // TM power reduction of 0.00136 is included in normalization as in Collin paper.
   // coupling*coupling is the power fraction plotted in the Locust paper.

        double coupling =   705.7 * 2./LMCConst::Pi() * 4./(2.*LMCConst::Pi()) / kc *
        		j1(kc*r);
      
        return coupling*coupling;
    }


    double FieldCalculator::GetCouplingFactorTE10(Kassiopeia::KSParticle& aFinalParticle)  // Phase 1
    {
        double dim1_wr42 = fInterface->fField->GetDimX(); // a in m
        double x = aFinalParticle.GetPosition().GetX() + dim1_wr42/2.;
        double coupling = 0.63*sin(LMCConst::Pi()*x/dim1_wr42);  // avg over cyclotron orbit.
        return coupling*coupling;
    }


    double FieldCalculator::GetTE10FieldAfterOneBounce(Kassiopeia::KSParticle& aFinalParticle)
    {
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTE10(aFinalParticle);
        double zvelocity = aFinalParticle.GetVelocity().GetZ();
        double zPosition = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocityTE10(aFinalParticle),2.),0.5);

        double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);
        double phi_shortTE10 = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(fabs(zPosition) + fInterface->fCENTER_TO_SHORT)/(GroupVelocity/fprime_short);  // phase of reflected field at position of electron.

        double FieldFromShort = cos(0.) + cos(phi_shortTE10);

	//        printf("field sum at z=%g is %f with zvelocity %g\n", zPosition, FieldFromShort, zvelocity); 
        //getchar();

        return FieldFromShort;  // Phase 1

    }



    double FieldCalculator::GetTM01FieldWithTerminator(Kassiopeia::KSParticle& aFinalParticle)
    {
        double tCyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocityTM01(aFinalParticle);
        double tVelocityZ = aFinalParticle.GetVelocity().GetZ();
        double tPositionZ = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/sqrt(1.0-pow(tVelocityZ / GetGroupVelocityTM01(aFinalParticle),2.) );
        double fprime_polarizer = tCyclotronFrequency*GammaZ*(1.-tVelocityZ/GroupVelocity);

        double phi_polarizerTM01 = 2.*LMCConst::Pi()*(2.*(fInterface->fCENTER_TO_ANTENNA-fabs(tPositionZ)))/(GroupVelocity/fprime_polarizer);
        double TM01FieldWithTerminator = cos(0.) + cos(phi_polarizerTM01);
        //printf("TM01FieldWithTerminator is %f\n", TM01FieldWithTerminator);
        //getchar();
        return TM01FieldWithTerminator;


    }


    double FieldCalculator::GetDampingFactorPhase1(Kassiopeia::KSParticle& aFinalParticle)
    {
      double TE10FieldFromShort = 0.;
	  TE10FieldFromShort = GetTE10FieldAfterOneBounce(aFinalParticle);
      double A10squ = GetCouplingFactorTE10(aFinalParticle);
      double DampingFactorTE10 = 1. - A10squ + A10squ*TE10FieldFromShort*TE10FieldFromShort;  // = P'/P
      return DampingFactorTE10;
    }



    double FieldCalculator::GetDampingFactorPhase2(Kassiopeia::KSParticle& aFinalParticle)
    {
      double TM01FieldWithTerminator = 0.;
      TM01FieldWithTerminator = GetTM01FieldWithTerminator(aFinalParticle);
      double A01squ = GetCouplingFactorTM01(aFinalParticle);
      double DampingFactorTM01 = 1. - A01squ + A01squ*TM01FieldWithTerminator*TM01FieldWithTerminator;  // = P'/P
      return DampingFactorTM01;
    }

    double FieldCalculator::GetDampingFactorCavity(Kassiopeia::KSParticle& aFinalParticle)
    {
    	return GetDampingFactorPhase1(aFinalParticle);
    }


    double FieldCalculator::GetCavityFIRSample(std::vector<double> tKassParticleXP, bool BypassTF)
    {
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double orbitPhase = tKassParticleXP[6];  // radians
    	double cycFrequency = tKassParticleXP[7];  // rad/s
    	double vMag = pow(tVx*tVx + tVy*tVy,0.5);
    	double convolution = 0.0;

		// populate FIR filter with frequency for just this sample interval:
		for (int i=0; i < fInterface->nFilterBinsRequired; i++)
		{
			fInterface->FIRfrequencyBuffer[0].push_back(cycFrequency);  // rad/s
			fInterface->FIRfrequencyBuffer[0].pop_front();
		}

		std::deque<double>::iterator it = fInterface->FIRfrequencyBuffer[0].begin();
		while (it != fInterface->FIRfrequencyBuffer[0].end())
		{
			orbitPhase += (*it)*fInterface->dtFilter;

			if (*it != 0.)
			{
				fInterface->ElementFIRBuffer[0].push_back(cos(orbitPhase));
			}
			else
			{
				fInterface->ElementFIRBuffer[0].push_back(0.);
			}
			fInterface->ElementFIRBuffer[0].pop_front();

			*it++;
		}

		if ( !BypassTF )
		{
			convolution = fInterface->fTFReceiverHandler.ConvolveWithFIRFilter(fInterface->ElementFIRBuffer[0]);
		}
		else
		{
			convolution = 1.0;
		}

		return LMCConst::Q()*vMag*convolution;

    }


    std::vector<double> FieldCalculator::GetCavityNormalizedModeField(int l, int m, int n, std::vector<double> tLocation, bool TE, bool Electric)
      {
      	double tR = tLocation[0];
      	double tZ = tLocation[2];
      	std::vector<double> tField;
      	double normFactor = 0.;

      	if (TE)
      	{
      		if (Electric)  // Get the electric field, usually at the electron.
      		{
      			tField = fInterface->fField->TE_E(l,m,n,tR,0.,tZ,1);
      		}
      		else  // Get the magnetic field, nominally at e.g. a readout probe location.
      		{
      			tField = fInterface->fField->TE_H(l,m,n,tR,0.,tZ,1);
      		}
      		normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];
      	}
      	else
      	{
      		if (Electric)
      		{
      			tField = fInterface->fField->TM_E(l,m,n,tR,0.,tZ,1);
      		}
      		else
      		{
      			tField = fInterface->fField->TM_H(l,m,n,tR,0.,tZ,1);
      		}
      		normFactor = fInterface->fField->GetNormFactorsTM()[l][m][n];
      	}


  		auto it = tField.begin();

  		while (it != tField.end())
  		{
  			if (!isnan(*it))
  				(*it) *= normFactor;
  			*it++;
  		}

      	return tField;  // return normalized field.
      }



    std::vector<double> FieldCalculator::GetWaveguideNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP)
     {
    	// The l index is inert in the waveguide.
     	double tX = tKassParticleXP[0] * cos(tKassParticleXP[1]);
     	double tY = tKassParticleXP[0] * sin(tKassParticleXP[1]);
     	double fcyc = tKassParticleXP[7];
     	std::vector<double> tTE_E_electron = fInterface->fField->TE_E(m,n,tX,tY,fcyc);
 		double normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];

 		auto it = tTE_E_electron.begin();
 		while (it != tTE_E_electron.end())
 		{
 			if (!isnan(*it))
 				(*it) *= normFactor;
 			*it++;
 		}
     	return tTE_E_electron;  // return normalized field.
     }


    double FieldCalculator::GetCavityDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, bool IntermediateFile)
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


    double FieldCalculator::GetWaveguideDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized, bool IntermediateFile)
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







}

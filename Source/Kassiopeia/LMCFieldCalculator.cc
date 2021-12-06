
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
        double CutOffFrequency = SpeedOfLight * LMCConst::Pi() / 10.668e-3; // a in m
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
        double dim1_wr42 = 10.668e-3; // a in m
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

    	std::vector<double> tTE_E_normalized = GetCavityNormalizedModeField(0,1,1,aFinalParticle); // lmn 011 for now.
    	fInterface->dotProductFactor = GetCavityDotProductFactor(aFinalParticle, tTE_E_normalized);  // unit velocity \dot unit theta
//    	fInterface->dotProductFactor = 0.5;  // TO-DO:  Check dotProductFactor - should it be an instantaneous dot product?
    	fInterface->modeAmplitude = tTE_E_normalized.back();  // absolute E_theta at electron

    	fInterface->CavityFIRSample = GetCavityFIRSample(aFinalParticle, fInterface->nFilterBinsRequired, fInterface->dtFilter);

    	double DampingFactorCavity = 1.0; // TO-DO:  Insert power correction factor here.
    	return DampingFactorCavity;
    }

    double FieldCalculator::GetCavityFIRSample(Kassiopeia::KSParticle& aFinalParticle, int nFilterBinsRequired, double dtFilter)
    {
    	double tVx = aFinalParticle.GetVelocity().X();
    	double tVy = aFinalParticle.GetVelocity().Y();
    	double orbitPhase = calcOrbitPhase(tVx, tVy);  // radians
    	double fieldFrequency = 2.*LMCConst::Pi() * aFinalParticle.GetCyclotronFrequency();  // rad/s
    	double vMag = pow(tVx*tVx + tVy*tVy,0.5);
    	double convolution = 0.0;

		// populate FIR filter with frequency for just this sample interval:
		for (int i=0; i < nFilterBinsRequired; i++)
		{
			fInterface->FIRfrequencyBuffer[0].push_back(fieldFrequency);  // rad/s
			fInterface->FIRfrequencyBuffer[0].pop_front();
		}

		// populate entire FIR filter with current, using frequencies from recent previous samples:
		std::deque<double>::iterator it = fInterface->FIRfrequencyBuffer[0].begin();
		while (it != fInterface->FIRfrequencyBuffer[0].end())
		{
			orbitPhase += (*it)*dtFilter;
			if (*it != 0.)
			{
				fInterface->ElementFIRBuffer[0].push_back(LMCConst::Q()*vMag*cos(orbitPhase));
			}
			else
			{
				fInterface->ElementFIRBuffer[0].push_back(0.);
			}
			fInterface->ElementFIRBuffer[0].pop_front();
			*it++;
		}

		convolution=fInterface->fTFReceiverHandler.ConvolveWithFIRFilter(fInterface->ElementFIRBuffer[0]);

		// Make copies of the buffers to keep the contents arranged sequentially in Locust:
		fInterface->ElementFIRBufferCopy[0] = fInterface->ElementFIRBuffer[0];
		fInterface->FIRfrequencyBufferCopy[0] = fInterface->FIRfrequencyBuffer[0];

		return convolution;

    }


    double FieldCalculator::quadrantOrbitCorrection(double phase, double vx)
    {
    	double phaseCorrection = 0.;
    	if (((phase < 0.)&&(vx < 0.)) || ((phase > 0.)&&(vx > 0.)))
    		phaseCorrection = LMCConst::Pi();

    	return phaseCorrection;
    }


    double FieldCalculator::quadrantPositionCorrection(double phase, double x)
    {
    	double phaseCorrection = 0.;
    	if (((phase < 0.)&&(x < 0.)) || ((phase > 0.)&&(x < 0.)))
    		phaseCorrection = LMCConst::Pi();

    	return phaseCorrection;
    }


    double FieldCalculator::calcOrbitPhase(double vx, double vy)
    {
    	double phase = 0.;
    	if (fabs(vy) > 0.)
    		phase = atan(-vx/vy);
    	phase += quadrantOrbitCorrection(phase, vx);
//    	printf("phase is %g\n", phase*180./LMCConst::Pi()); getchar();
    	return phase;
    }

    double FieldCalculator::calcTheta(double x, double y)
    {
    	double phase = 0.;
    	if (fabs(x) > 0.)
    		phase = atan(y/x);
    	phase += quadrantPositionCorrection(phase, x);
    	return phase;
    }



    double FieldCalculator::GetCavityDotProductFactor(Kassiopeia::KSParticle& aFinalParticle, std::vector<double> aTE_E_normalized)
    {
    	double tX = aFinalParticle.GetPosition().X();
    	double tY = aFinalParticle.GetPosition().Y();
    	double tThetaParticle = calcTheta(tX, tY);

    	double tEtheta = aTE_E_normalized.back();
    	double tEx = -sin(tThetaParticle) * tEtheta;
    	double tEy = cos(tThetaParticle) * tEtheta;
    	double tEmag = tEtheta;
    	double tVx = aFinalParticle.GetVelocity().X();
    	double tVy = aFinalParticle.GetVelocity().Y();
    	double tVmag = pow(tVx*tVx + tVy*tVy, 0.5);
    	return fabs(tEx*tVx + tEy*tVy)/tEmag/tVmag;  // fabs ( unit J \cdot unit E )
    }



    std::vector<double> FieldCalculator::GetCavityNormalizedModeField(int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle)
    {
    	double tX = aFinalParticle.GetPosition().X();
    	double tY = aFinalParticle.GetPosition().Y();
    	double tZ = aFinalParticle.GetPosition().Z();
    	double tR = sqrt(tX*tX + tY*tY);
    	double fcyc = aFinalParticle.GetCyclotronFrequency();  // TO_DO Is this in radians or Hz?
    	std::vector<double> tTE_E_electron = fInterface->fField->TE_E(l,m,n,tR,0.,tZ, fcyc);
		double normFactor = fInterface->fField->GetNormFactorsTE()[l][m][n];  // select mode 0,1,1

		auto it = tTE_E_electron.begin();
		while (it != tTE_E_electron.end())
		{
			if (!isnan(*it))
				(*it) *= normFactor;
			*it++;
		}
    	return tTE_E_electron;  // return normalized field.
    }



}


#include "LMCFieldCalculator.hh"
#include "logger.hh"
#include "KSParticleFactory.h"
#include <algorithm>

#include "KSInteractionsMessage.h"
#include <limits>

using std::numeric_limits;

using namespace Kassiopeia;
namespace locust
{

    LOGGER( lmclog, "FieldCalculator" );

    FieldCalculator::FieldCalculator() :
    		fNFilterBinsRequired(0),
			fTFReceiverHandler( NULL ),
			fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }
    FieldCalculator::FieldCalculator( const FieldCalculator& aCopy ) :
    		fNFilterBinsRequired(0),
			fTFReceiverHandler( NULL ),
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

    bool FieldCalculator::ConfigureByInterface()
    {
    	if (fInterface->fConfigureKass)
    	{
    	    const scarab::param_node* aParam = fInterface->fConfigureKass->GetParameters();
    	    if (!this->Configure( *aParam ))
    	    {
    		    LERROR(lmclog,"Error configuring FieldInterface class");
    		    exit(-1);
    	    }
    	}
        return true;
    }

    bool FieldCalculator::Configure( const scarab::param_node& aParam )
     {

    	fTFReceiverHandler = new TFReceiverHandler;
    	if(!fTFReceiverHandler->Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver FIRHandler class");
    	}

        if( aParam.has( "tf-receiver-filename" ) )
        {
            if (!fTFReceiverHandler->ReadHFSSFile())  // Read external file
            {
            	LERROR(lmclog,"FIR has not been generated.");
            	exit(-1);
            }
        }
        else // Generate analytic response function
        {
        	if ((aParam.has( "equivalent-circuit" ) ) && (aParam["equivalent-circuit"]().as_bool()))
        	{
        		fAnalyticResponseFunction = new EquivalentCircuit();
        		if ( !fAnalyticResponseFunction->Configure(aParam) )
        		{
        			LWARN(lmclog,"EquivalentCircuit was not configured.");
        			return false;
        		}
        		else
        		{
        			if (!fTFReceiverHandler->ConvertAnalyticTFtoFIR(fAnalyticResponseFunction->GetInitialFreq(),fAnalyticResponseFunction->GetTFarray()))
        			{
        				LWARN(lmclog,"TF->FIR was not generated correctly.");
        				return false;
        			}
        		}
        	}
        	else
        	{
        		fAnalyticResponseFunction = new DampedHarmonicOscillator();
        		if ( !fAnalyticResponseFunction->Configure(aParam) )
        		{
        			LWARN(lmclog,"DampedHarmonicOscillator was not configured.");
        			return false;
        		}
        		if (!fTFReceiverHandler->ConvertAnalyticGFtoFIR(0,1,1,fAnalyticResponseFunction->GetGFarray(0,1,1)))
        		{
        			LWARN(lmclog,"GF->FIR was not generated.");
        			return false;
        		}
        	}
        } // aParam.has( "tf-receiver-filename" )

        SetFilterSize( fTFReceiverHandler->GetFilterSize());

        return true;
     }

    void FieldCalculator::SetNFilterBinsRequired( double dt )
    {
    	if (fTFReceiverHandler)
    	{
    	    fNFilterBinsRequired = 1 + (int)( (dt) / fTFReceiverHandler->GetFilterResolution());
	std::cout << "fNFilterBinsRequired set to " << fNFilterBinsRequired << std::endl;
    	}
    }

    int FieldCalculator::GetNFilterBinsRequired()
    {
    	return fNFilterBinsRequired;
    }

    void FieldCalculator::SetFilterSize( int aFilterSize )
    {
    	fFIRBuffer.resize( aFilterSize );
    	fFrequencyBuffer.resize( aFilterSize );
    }

    int FieldCalculator::GetFilterSize()
    {
    	return fFIRBuffer.size();
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

    double FieldCalculator::GetCouplingFactorCavity(int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle)
    {
    	double dimR = fInterface->fField->GetDimR(); // m
    	double x = aFinalParticle.GetPosition().GetX();
    	double y = aFinalParticle.GetPosition().GetY();
    	double z = aFinalParticle.GetPosition().GetZ();
    	double r = pow( x*x + y*y, 0.5 );
    	double norm = fInterface->fField->GetTE_E(l,m,n,dimR/2.,0.,0.,false).back(); // max value
    	double tAvgDotProductFactor = 0.63;  // TO-DO:  this should be calculated and not just overridden like this.
    	double coupling = tAvgDotProductFactor * fInterface->fField->GetTE_E(l,m,n,r,0.,z,false).back()/norm;
    	return coupling*coupling;
    }

    double FieldCalculator::GetFieldCavity(int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle)
    {
    	std::vector<double> tKassParticleXP = fInterface->fTransmitter->ExtractParticleXP(0., false, 0., false);
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double vMag = pow(tVx*tVx + tVy*tVy,0.5);
	std::cout << "About to convolve" << std::endl;
        std::pair<double,double> complexConvolution = GetCavityFIRSample(l,m,n,tKassParticleXP, 0);
	std::cout << "Convolved with elements: " << std::endl;
	std::cout << complexConvolution.first << " " << complexConvolution.second << std::endl;
        // The excitation amplitude A_\lambda should be calculated the same way here
        // as in the signal generator, but here it is normalized such that the
        // E-field magnitude peaks at the Hanneke factor when the cyclotron frequency
        // is on the mode resonance.

        // Convolution with LMCDampedHarmonicOscillator resonance peaks at 1.0*fHannekePowerFactor,
        // and GetCavityFIRSample returns that convolution scaled with vMag*Q.  Normalizing with
        // vMag*Q as below, we have a dhoMag that peaks at 1.0*fHannekePowerFactor:

        double dhoNorm = vMag * LMCConst::Q();
        double dhoMag = 0.;
        if (dhoNorm > 0.)
        {
        	dhoMag = complexConvolution.first / dhoNorm;
        }
        double dhoPhase = complexConvolution.second;

        // first term represents the new field driven by the electron.
        // second term represents the field driven by the electron previously.
        // "dhoMag" scales from ~0 (off-resonance) to DampedHarmonicOscillator::fHannekePowerFactor (on-resonance).
        // The electron should radiate maximally if on resonance.
        double fieldCavity = cos(0.) + dhoMag*cos(dhoPhase);
        return fieldCavity;
    }


    double FieldCalculator::GetDampingFactorCavity(int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle)
    {
    	double FieldFromCavity = GetFieldCavity(l,m,n,aFinalParticle);
    	double Asqu = GetCouplingFactorCavity(l,m,n,aFinalParticle);
    	double DampingFactorCavity = 1. - Asqu + Asqu*FieldFromCavity*FieldFromCavity;  // = P'/P
    	if ( fabs(DampingFactorCavity) > 0. )
    	{
    		return DampingFactorCavity;
    	}
    	else
    	{
    		return 1.0;
    	}
    }


    std::pair<double,double> FieldCalculator::GetCavityFIRSample(int l, int m, int n, std::vector<double> tKassParticleXP, bool BypassTF)
    {
    	double convolutionMag = 0.0;
    	double convolutionPhase = 0.0;
    	double tVx = tKassParticleXP[3];
    	double tVy = tKassParticleXP[4];
    	double vMag = pow(tVx*tVx + tVy*tVy,0.5);

    	if ( !BypassTF )
    	{
    		double orbitPhase = tKassParticleXP[6];  // radians
    		double cycFrequency = tKassParticleXP[7];  // rad/s
    		// populate FIR filter with frequency for just this sample interval:
		std::cout << "fNFilterBinsRequired: " << std::endl;
		std::cout << fNFilterBinsRequired << std::endl;
    		for (int i=0; i < fNFilterBinsRequired; i++)
    		{
    			fFrequencyBuffer.push_back(cycFrequency);  // rad/s
    			fFrequencyBuffer.pop_front();
    		}

    		std::deque<double>::iterator it = fFrequencyBuffer.begin();
		std::cout << "fFrequencyBuffer size: " << fFrequencyBuffer.size() << std::endl;
    		while (it != fFrequencyBuffer.end())
    		{
    			// TO-DO:  Consider:  Replace dtFilter with z(t)/vp.
    			orbitPhase += (*it)*fTFReceiverHandler->GetFilterResolution();

    			if (*it != 0.)
    			{
    				fFIRBuffer.push_back(cos(orbitPhase));
    			}
    			else
    			{
    				fFIRBuffer.push_back(0.);
    			}
    			fFIRBuffer.pop_front();

    			*it++;
    		}
		std::cout << "fFIRBuffer size: " << fFIRBuffer.size() << std::endl;
                std::pair<double,double> convolution = fTFReceiverHandler->ConvolveWithComplexFIRFilterArray(l,m,n,fFIRBuffer);
    		convolutionMag = convolution.first;
    		convolutionPhase = convolution.second;

    	}
    	else
    	{
    		convolutionMag = 1.0;
    		convolutionPhase = 0.;
    	}

    	return std::make_pair(convolutionMag*LMCConst::Q()*vMag, convolutionPhase);

    }





}

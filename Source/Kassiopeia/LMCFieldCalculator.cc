
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
        fNFilterBinsRequired( 0 ),
        fTFReceiverHandler( NULL ),
        fAnalyticResponseFunction( 0 ),
        fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }
    FieldCalculator::FieldCalculator( const FieldCalculator& aCopy ) :
        fNFilterBinsRequired( 0 ),
        fTFReceiverHandler( NULL ),
        fAnalyticResponseFunction( 0 ),
        fInterface( aCopy.fInterface )
    {
    }
    FieldCalculator* FieldCalculator::Clone() const
    {
        return new FieldCalculator( *this );
    }
    FieldCalculator::~FieldCalculator()
    {
        if (fTFReceiverHandler != NULL)
        {
            delete fTFReceiverHandler;
        }
        if (fAnalyticResponseFunction != NULL)
        {
            delete fAnalyticResponseFunction;
        }
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
        fModeSet = fInterface->fField->ModeSelect(fInterface->fbWaveguide, 0);

        fTFReceiverHandler = new TFReceiverHandler;
        if(!fTFReceiverHandler->Configure(aParam))
        {
            LERROR(lmclog,"Error configuring receiver FIRHandler class");
        }

        if( aParam.has( "tf-receiver-filename" ) )
        {
            if (!fTFReceiverHandler->ReadHFSSFile(1,0,1,0))  // Read external file
            {
                LERROR(lmclog,"FIR has not been generated.");
                exit(-1);
            }
        }
        else // Generate analytic response function
        {
            fAnalyticResponseFunction = new DampedHarmonicOscillator();
            if ( !fAnalyticResponseFunction->Configure(aParam) )
            {
                LWARN(lmclog,"DampedHarmonicOscillator was not configured.");
                return false;
            }
            if (!fTFReceiverHandler->ConvertAnalyticGFtoFIR(fModeSet, fAnalyticResponseFunction->GetGFarray(fModeSet)))
            {
                LWARN(lmclog,"GF->FIR was not generated.");
                return false;
            }
        } // aParam.has( "tf-receiver-filename" )

        // Size the electron information buffers to be similar to the Green's functions:
        SetFilterSize( fTFReceiverHandler->GetFilterSizeArray(fModeSet[0][0],fModeSet[0][1],fModeSet[0][2],fModeSet[0][3]));
        return true;
    }

    void FieldCalculator::SetNFilterBinsRequired( double dt )
    {
        if (fTFReceiverHandler)
        {
            fNFilterBinsRequired = 1 + (int)( (dt) / fTFReceiverHandler->GetFilterResolutionArray(fModeSet[0][0],fModeSet[0][1],fModeSet[0][2],fModeSet[0][3]));
        }
    }

    int FieldCalculator::GetNFilterBinsRequired()
    {
        return fNFilterBinsRequired;
    }

    void FieldCalculator::SetFilterSize( int aFilterSize )
    {
        // These contain histories of the electron's orbit phase and cyclotron frequency:
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
        double phi_shortTE10 = LMCConst::Pi()/2. + 2.*LMCConst::Pi()*(fabs(zPosition) + fInterface->fField->GetCenterToShort())/(GroupVelocity/fprime_short);  // phase of reflected field at position of electron.

        double FieldFromShort = cos(0.) + cos(phi_shortTE10);

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

        double phi_polarizerTM01 = 2.*LMCConst::Pi()*(2.*(fInterface->fField->GetCenterToAntenna()-fabs(tPositionZ)))/(GroupVelocity/fprime_polarizer);
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

    double FieldCalculator::GetCouplingFactorTXlmnCavity(int l, int m, int n, bool bTE, Kassiopeia::KSParticle& aFinalParticle)
    {
        double tAvgDotProductFactor = fInterface->fField->GetAvgDotProductFactor()[l][m][n];
        double norm = 0.;
        double coupling = 0.;
        double dimR = fInterface->fField->GetDimR(); // m
        double dimZ = fInterface->fField->GetDimL(); // m
        double x = aFinalParticle.GetPosition().GetX();
        double y = aFinalParticle.GetPosition().GetY();
        double z = aFinalParticle.GetPosition().GetZ();
        double r = pow( x*x + y*y, 0.5 );
        double tZ = dimZ/2. - dimZ/2./n;  // Normalize in an antinode along z

        if (bTE)
        {
    	    norm = fInterface->fField->GetTE_E(l,m,n,dimR/2.,0.,tZ,false).back(); // max value, TO-DO:  make more general?
    	    coupling = tAvgDotProductFactor * fInterface->fField->GetTE_E(l,m,n,r,0.,z,false).back()/norm;
        }
        else
        {
            norm = fInterface->fField->GetTM_E(l,m,n,dimR/2.,0.,tZ,false).back(); // TO-DO:  decide coordinates
            coupling = tAvgDotProductFactor * fInterface->fField->GetTM_E(l,m,n,r,0.,z,false).back()/norm;
        }

        return coupling*coupling;
    }

    double FieldCalculator::GetTXlmnFieldCavity(int l, int m, int n, bool bTE, Kassiopeia::KSParticle& aFinalParticle)
    {

        // l, m, & n are needed for selecting the resonant frequency and Q.  (Still TO-DO).

        double tVx = aFinalParticle.GetVelocity().X();
        double tVy = aFinalParticle.GetVelocity().Y();

        std::vector<double> tKassParticleXP;
        tKassParticleXP.push_back(aFinalParticle.GetPosition().X());
        tKassParticleXP.push_back(aFinalParticle.GetPosition().Y());
        tKassParticleXP.push_back(aFinalParticle.GetPosition().Z());
        tKassParticleXP.push_back(aFinalParticle.GetVelocity().X());
        tKassParticleXP.push_back(aFinalParticle.GetVelocity().Y());
        tKassParticleXP.push_back(aFinalParticle.GetVelocity().Z());
        tKassParticleXP.push_back(calcOrbitPhase(tVx, tVy));
        tKassParticleXP.push_back(aFinalParticle.GetCyclotronFrequency() * 2. * LMCConst::Pi());
        tKassParticleXP.push_back(aFinalParticle.GetTime());

        double vMag = pow(tVx*tVx + tVy*tVy,0.5);

        std::pair<double,double> complexConvolution = GetCavityFIRSample(bTE, l, m, n, tKassParticleXP, 0);

        // The excitation amplitude A_\lambda should be calculated the same way here
        // as in the signal generator.

        // Convolution with LMCDampedHarmonicOscillator resonance peaks at 1.0,
        // and GetCavityFIRSample returns that convolution scaled with vMag*Q.  Normalizing with
        // vMag*Q as below, we have a dhoMag that peaks at 1.0:

        double dhoNorm = vMag * LMCConst::Q();
        double dhoMag = 0.;
        if (dhoNorm > 0.)
        {
            dhoMag = complexConvolution.first / dhoNorm;
        }

        // For compatibility with the develop branch, keep dhoPhase in the 1st and 4th
        // quadrants with atan(tan(dhoPhase)).  This offsets the pre-existing field phase from
        // the drive phase.  This may need to be revisited by removing the atan(tan()) for a
        // closer match between the pre-existing field phase and the drive phase.  The track
        // slope is sensitive to this phase.
        double dhoPhase = atan(tan(complexConvolution.second));

        // first term represents the new field driven by the electron.
        // second term represents the field driven by the electron previously.
        // "dhoMag" scales from ~0 (off-resonance) to DampedHarmonicOscillator::fHannekePowerFactor (on-resonance).
        // The electron should radiate maximally if on resonance.
        double fieldCavity = cos(0.) + dhoMag*cos(dhoPhase);

        return fieldCavity;
    }


    double FieldCalculator::GetDampingFactorCavity(Kassiopeia::KSParticle& aFinalParticle)
    {
        double DampingFactorCavity = 0.;
        for (int mu=0; mu<fModeSet.size(); mu++)
        {
            bool bTE = fModeSet[mu][0];
            int l = fModeSet[mu][1];
            int m = fModeSet[mu][2];
            int n = fModeSet[mu][3];

            double TXlmnFieldFromCavity = GetTXlmnFieldCavity(l,m,n,bTE,aFinalParticle);
            double Almnsqu = GetCouplingFactorTXlmnCavity(l,m,n,bTE,aFinalParticle);
            double DampingFactorTXlmnCavity = 1. - Almnsqu + Almnsqu*TXlmnFieldFromCavity*TXlmnFieldFromCavity;  // = (P'/P)_{lmn}
            DampingFactorCavity += DampingFactorTXlmnCavity - 1.; // (P'/P)_{lmn} - 1
        }
    	if (fabs(DampingFactorCavity) > 0.)
    		return DampingFactorCavity + 1.0;
    	else
    		return 1.0;  // No feedback
    }

    std::pair<double,double> FieldCalculator::GetCavityFIRSample(int bTE, int l, int m, int n, std::vector<double> tKassParticleXP, bool BypassTF)
    {
        double convolutionMag = 0.0;
        double convolutionPhase = 0.0;
        double tVx = tKassParticleXP[3];
        double tVy = tKassParticleXP[4];
        double vMag = pow(tVx*tVx + tVy*tVy,0.5);
        double orbitPhase = tKassParticleXP[6];  // radians
        double cycFrequency = tKassParticleXP[7];
        double tTime = tKassParticleXP[9];
        double amplitude = 0.;
        if ( fInterface->fField->InVolume(tKassParticleXP))
        {
            amplitude = 1.;
        }


        if ( !BypassTF )
        {
            // populate FIR filter with frequency for just this sample interval:
            for (int i=0; i < fNFilterBinsRequired; i++)
            {
                fFrequencyBuffer.push_back(cycFrequency);  // rad/s
                fFrequencyBuffer.pop_front();
            }

            std::deque<double>::iterator it = fFrequencyBuffer.begin();
            while (it != fFrequencyBuffer.end())
            {
                orbitPhase += (*it)*fTFReceiverHandler->GetFilterResolutionArray(bTE, l, m, n);

                if (*it != 0.)
                {
                    fFIRBuffer.push_back(amplitude * cos(orbitPhase));
                }
                else
                {
                    fFIRBuffer.push_back(0.);
                }
                fFIRBuffer.pop_front();

                *it++;
            }

            std::pair<double,double> convolution = fTFReceiverHandler->ConvolveWithComplexFIRFilterArray(bTE, l, m, n, fFIRBuffer);

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

    double FieldCalculator::calcOrbitPhase(double vx, double vy)
    {
        double phase = 0.;
        if ((fabs(vy) > 0.))
        {
            phase = atan(-vx/vy);
        }

        phase += quadrantOrbitCorrection(phase, vx);
        return phase;
    }

    double FieldCalculator::quadrantOrbitCorrection(double phase, double vx)
    {
        double phaseCorrection = 0.;
        if (((phase < 0.)&&(vx < 0.)) || ((phase > 0.)&&(vx > 0.)))
            phaseCorrection = LMCConst::Pi();

        return phaseCorrection;
    }



}

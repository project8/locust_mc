
#include "LMCPowerNormFieldCalculator.hh"
#include "logger.hh"
#include "KSParticleFactory.h"
#include <algorithm>

#include "KSInteractionsMessage.h"
#include <limits>

using std::numeric_limits;
using katrin::KThreeVector;


using namespace Kassiopeia;
namespace locust
{

    LOGGER( lmclog, "PowerNormFieldCalculator" );

    PowerNormFieldCalculator::PowerNormFieldCalculator() :
        fNFilterBinsRequired( 0 ),
        fTFReceiverHandler( NULL ),
        fAnalyticResponseFunction( 0 ),
        fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }
    PowerNormFieldCalculator::PowerNormFieldCalculator( const PowerNormFieldCalculator& aCopy ) :
        fNFilterBinsRequired( 0 ),
        fTFReceiverHandler( NULL ),
        fAnalyticResponseFunction( 0 ),
        fInterface( aCopy.fInterface )
    {
    }
    PowerNormFieldCalculator* PowerNormFieldCalculator::Clone() const
    {
        return new PowerNormFieldCalculator( *this );
    }
    PowerNormFieldCalculator::~PowerNormFieldCalculator()
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

    bool PowerNormFieldCalculator::Configure( const scarab::param_node& aParam )
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
        SetFilterSize( 100 );
        // Set the initial time to zero. TODO: this should be set to the time the particle is emitted
        fTime = 0.0;

        return true;
    }

    void PowerNormFieldCalculator::SetNFilterBinsRequired( double dt )
    {
        if (fTFReceiverHandler)
        {
            fNFilterBinsRequired = 1 + (int)( (dt) / fTFReceiverHandler->GetFilterResolutionArray(fModeSet[0][0],fModeSet[0][1],fModeSet[0][2],fModeSet[0][3]));
        }
    }

    double PowerNormFieldCalculator::GetDampingFactorCavity(Kassiopeia::KSParticle& aFinalParticle)
    {
        double deltaE = 0.;
        for (int mu=0; mu<fModeSet.size(); mu++)
        {
            bool bTE = fModeSet[mu][0];
            int l = fModeSet[mu][1];
            int m = fModeSet[mu][2];
            int n = fModeSet[mu][3];

            SetCavityFIRSample(bTE, l, m, n, aFinalParticle, 0);

            const std::vector<std::array<double, 2>>& tEfield = fTFReceiverHandler->GetEfield()[bTE][l][m][n];

            // LPROG("B field at set: " << fTFReceiverHandler->GetBfield()[bTE][l][m][n].back()[0]);

            int timestep = 0;

            double dt = fTFReceiverHandler->GetFilterResolutionArray(bTE, l, m, n);

            for (auto it = fJdotEBuffer.begin();it!=fJdotEBuffer.end(); ++it)
            {
                deltaE += (*it) * tEfield[timestep][0] * dt;
                timestep++;
            }
        }
        return deltaE / 1e7;
    }

    void PowerNormFieldCalculator::SetCavityFIRSample(int bTE, int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle, bool BypassTF)
    {
        KThreeVector tPosition = aFinalParticle.GetPosition();
        KThreeVector tVelocity = aFinalParticle.GetVelocity();
        KThreeVector tGuidingCenterPosition = aFinalParticle.GetGuidingCenterPosition();
        KThreeVector tMagneticField = aFinalParticle.GetMagneticField();
        double tCyclotronFrequency = - aFinalParticle.GetCyclotronFrequency(); // - sign for electrons only

        KThreeVector tVelocityParallel = tVelocity.Dot(tMagneticField.Unit()) * tMagneticField.Unit();
        KThreeVector vPerp = tVelocity - tVelocityParallel;

        double tCyclotronRadius = - vPerp.Magnitude() / tCyclotronFrequency / 2 / LMCConst::Pi();

        KThreeVector tBeta = vPerp.Unit();
        KThreeVector tAlpha = tPosition - tGuidingCenterPosition;
        tAlpha = tAlpha.Unit();

        double tTime = aFinalParticle.GetTime();

        // Number of time steps between tTime and fTime
        double dt = fTFReceiverHandler->GetFilterResolutionArray(bTE, l, m, n);
        int Nsteps = std::round((tTime - fTime) / dt);

        fJdotEBuffer.clear();
        // populate FIR filter with frequency for just this sample interval:
        // LPROG("tTime: " << tTime << "tTime - Nsteps*dt: " << tTime - Nsteps*dt << " dt: " << dt);
        // LPROG("Cyclotron frequency: " << tCyclotronFrequency);
        // if (tTime > 4.e-10) exit(-1);
        for (int i=0; i < Nsteps; i++)
        {
            // Test
            // fJdotEBuffer.push_back(LMCConst::Q() * LMCConst::C() * sin(2 * LMCConst::Pi() * tCyclotronFrequency * (fTime + dt*i)) * 1000);
            // fJdotEBuffer.push_back(LMCConst::Q() * LMCConst::C() * sin(2 * LMCConst::Pi() * 25.9e9 * (fTime + dt*i)) * 100);
            // LPROG("JdotEBuffer[" << i << "]: " << fJdotEBuffer[i]);

            // Calculate interpolated velocity
            KThreeVector interpolatedVelocity = InterpolateVelocity(-static_cast<double>(Nsteps-i)*dt, tCyclotronFrequency, tCyclotronRadius, tVelocityParallel, tMagneticField, tAlpha, tBeta);

            // Calculate interpolated position
            KThreeVector interpolatedPosition = InterpolatePosition(-static_cast<double>(Nsteps-i)*dt, tCyclotronFrequency, tCyclotronRadius, tVelocityParallel, tMagneticField, tGuidingCenterPosition, tAlpha, tBeta);

            // Create and populate particle state vector
            double thisR = pow(interpolatedPosition.X()*interpolatedPosition.X() + interpolatedPosition.Y()*interpolatedPosition.Y(), 0.5);
            // double thisTheta = calcTheta(interpolatedPosition.X(), interpolatedPosition.Y());
            double thisTheta = atan2(interpolatedPosition.X(), interpolatedPosition.Y());
            std::vector<double> tKassParticleXP;
            tKassParticleXP.push_back(thisR);
            tKassParticleXP.push_back(thisTheta);
            tKassParticleXP.push_back(interpolatedPosition.Z());
            tKassParticleXP.push_back(interpolatedVelocity.X());
            tKassParticleXP.push_back(interpolatedVelocity.Y());
            tKassParticleXP.push_back(interpolatedVelocity.Z());

            if ( !BypassTF && fInterface->fField->InVolume(tKassParticleXP))
            {
                // Calculate JdotE
                std::vector<double> tE_normalized = fInterface->fField->GetNormalizedModeField(l,m,n,tKassParticleXP,1,bTE);
                double tAvgDotProductFactor = fInterface->fField->CalculateDotProductFactor(l, m, n, tKassParticleXP, tE_normalized, 0);

                // Add to buffer
                fJdotEBuffer.push_back(tAvgDotProductFactor);
            }
            // if (i==0)
            // {
            //     LPROG("tKassParticleXP: ");
            //     for (size_t j = 0; j < tKassParticleXP.size(); ++j) {
            //         LPROG("  Element " << j << std::setprecision(16) << ": " << tKassParticleXP[j]);
            //     }
            // }
            tKassParticleXP.clear();

            // LPROG("JdotEBuffer[" << i << "]: " << fJdotEBuffer[i]);
        }
        // LPROG("Nsteps: " << Nsteps);
        // LPROG("Kass final time: " << std::setprecision(16) << tTime << " Sim time: " << tTime - Nsteps*dt);
        // LPROG("JdotEBuffer size: " << fJdotEBuffer.size());
        // LPROG("First JdotE: " << fJdotEBuffer[0]);

        fTFReceiverHandler->ComputeFields(bTE, l, m, n, fJdotEBuffer, fTime);

        // t_off = ; // insure that there is a time step dt between function calls
        fTime = tTime; // update the time
    }

    std::pair<double,double> PowerNormFieldCalculator::GetCavityFIRSample(int bTE, int l, int m, int n)
    {
        double Ereal = fTFReceiverHandler->GetEfield()[bTE][l][m][n].back()[0];
        double Breal = fTFReceiverHandler->GetBfield()[bTE][l][m][n].back()[0];

        return std::make_pair(Breal, Ereal);

    }

    KThreeVector PowerNormFieldCalculator::InterpolateVelocity(double dt, double tCyclotronFrequency, double tCyclotronRadius, KThreeVector tVelocityParallel, KThreeVector fMagneticField, KThreeVector tAlpha, KThreeVector tBeta)
    {
        KThreeVector tNewVelocity = tVelocityParallel + 2 * LMCConst::Pi() * tCyclotronRadius * tCyclotronFrequency *( - sin(2 * LMCConst::Pi() * tCyclotronFrequency * dt) * tAlpha + cos( 2 * LMCConst::Pi() * tCyclotronFrequency * dt ) * tBeta );

        return tNewVelocity;
    }

    KThreeVector PowerNormFieldCalculator::InterpolatePosition(double dt, double tCyclotronFrequency, double tCyclotronRadius, KThreeVector tVelocityParallel, KThreeVector tMagneticField, KThreeVector tGuidingCenterPosition, KThreeVector tAlpha, KThreeVector tBeta)
    {
        KThreeVector tNewPosition = tGuidingCenterPosition + tVelocityParallel * dt + tCyclotronRadius * ( cos( 2 * LMCConst::Pi() * tCyclotronFrequency * dt) * tAlpha + sin( 2 * LMCConst::Pi() * tCyclotronFrequency * dt) * tBeta);

        return tNewPosition;
    }




}

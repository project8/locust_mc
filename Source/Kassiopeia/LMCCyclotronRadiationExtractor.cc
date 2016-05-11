/*
 * LMCCyclotronRadiationExtractor.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCCyclotronRadiationExtractor.hh"
#include "LMCGlobalsDeclaration.hh"

using namespace Kassiopeia;
namespace locust
{

    CyclotronRadiationExtractor::CyclotronRadiationExtractor()
    {
    }

    CyclotronRadiationExtractor::CyclotronRadiationExtractor( const CyclotronRadiationExtractor& aOrig )
    {
    }

    CyclotronRadiationExtractor::~CyclotronRadiationExtractor()
    {
    }

    CyclotronRadiationExtractor* CyclotronRadiationExtractor::Clone() const
    {
        return new CyclotronRadiationExtractor( *this );
    }


    bool CyclotronRadiationExtractor::ExecutePreStepModification( KSParticle& anInitialParticle, KSParticleQueue& aQueue )
    {

    	return true;

    }

    bool CyclotronRadiationExtractor::ExecutePostStepModifcation( KSParticle& anInitialParticle, KSParticle& aFinalParticle, KSParticleQueue& aQueue )
    {


    	t_poststep = aFinalParticle.GetTime();

        if (t_poststep - t_old > 5.e-9)
            {
        	std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
            tLock.lock();
            Z = aFinalParticle.GetPosition().Z();
            X = aFinalParticle.GetPosition().X();
            de = aFinalParticle.GetKineticEnergy_eV() - anInitialParticle.GetKineticEnergy_eV();
            dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();
            xvelocity = aFinalParticle.GetVelocity().GetX();
            yvelocity = aFinalParticle.GetVelocity().GetY();
            zvelocity = aFinalParticle.GetVelocity().GetZ();

            GammaZ = 1.0/pow(1.0-pow(zvelocity/2.99792e8,2.),0.5);  // fix speed of light.

            fcyc = aFinalParticle.GetCyclotronFrequency();
            LarmorPower = -de/dt*1.602677e-19;
            tLock.unlock();

 //           printf("z position is %g and dt is %g\n", Z, dt);

            printf("Kassiopeia is about to send out a tick\n");
            fDigitizerCondition.notify_one();  // notify Locust after writing.

//                getchar();

             t_old = t_poststep;
             printf("de is %g and dt is %g and LarmorPower is %g\n", de, dt, LarmorPower);

          	 printf("Kassiopeia says:  tick has happened; continuous time is %g and zvelocity is %f\n", t_poststep, zvelocity);
          	 printf("Kassiopeia says:  fcyc is %g\n", fcyc);
          	 printf("  initial particle momentum: %g\n", anInitialParticle.GetMomentum().Z());
             printf("Mass is %g\n", anInitialParticle.GetMass());
             printf("k.e. = %g eV\n", anInitialParticle.GetKineticEnergy_eV());
             printf("Lorentz factor is %f\n", anInitialParticle.GetLorentzFactor());
             printf("1/(sqrt(1-v^2/c^2) is %f\n", 1.0/pow(1.0-pow(anInitialParticle.GetSpeed()/2.99792e8,2.),0.5));
             printf("FinalParticle().IsActive() is %d\n", aFinalParticle.IsActive());
          	}



        return true;
    }


    void CyclotronRadiationExtractor::InitializeComponent()
    {
    }

    void CyclotronRadiationExtractor::DeinitializeComponent()
    {
    }

    void CyclotronRadiationExtractor::PullDeupdateComponent()
    {
    }
    void CyclotronRadiationExtractor::PushDeupdateComponent()
    {
    }


} /* namespace locust */

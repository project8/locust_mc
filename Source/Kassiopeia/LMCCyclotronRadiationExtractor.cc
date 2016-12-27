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

        double aTimeStep = 0.;
    	t_poststep = aFinalParticle.GetTime();

        if (t_poststep - t_old > 5.e-10)
            {
        	std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
            tLock.lock();

        	KSTrajInterpolatorHermite anInterpolator;
        	KSTrajAdiabaticParticle anIntermediateAdiabaticParticle;
        	KSTrajAdiabaticParticle anInitialAdiabaticParticle;
        	KSTrajAdiabaticParticle aFinalAdiabaticParticle;
        	KSTrajAdiabaticIntegrator* anIntegrator;
        	KSTrajAdiabaticDifferentiator* aDifferentiator;


        	anInitialAdiabaticParticle.PullFrom(anInitialParticle);
        	aFinalAdiabaticParticle.PullFrom(aFinalParticle);
        	aTimeStep = t_old + 5.e-10 - anInitialParticle.GetTime();
//        	printf("aTimeStep is %g\n", aTimeStep); getchar();

            anInterpolator.GetInterpolate(0., *anIntegrator, *aDifferentiator, anInitialAdiabaticParticle, aFinalAdiabaticParticle, aTimeStep, anIntermediateAdiabaticParticle);


            Z = anIntermediateAdiabaticParticle.GetPosition().Z();
            X = anIntermediateAdiabaticParticle.GetPosition().X();
            de = anIntermediateAdiabaticParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy_eV();
            dt = anIntermediateAdiabaticParticle.GetTime() - anInitialParticle.GetTime();  // for step.  not digitizer.
            fcyc = anIntermediateAdiabaticParticle.GetCyclotronFrequency();
            t_old += 5.e-10;

/*
//                        Z = aFinalParticle.GetPosition().Z();
//                        X = aFinalParticle.GetPosition().X();
                        de = aFinalParticle.GetKineticEnergy_eV() - anInitialParticle.GetKineticEnergy_eV();
                        dt = t_poststep - t_old;
//                        fcyc = aFinalParticle.GetCyclotronFrequency();
                        t_old = t_poststep;
                        */


//printf("Zold was %g and Zinterpolated is %g\n", anInitialParticle.GetPosition().Z(), Z);
//            printf("dt is %g\n", dt);
//            getchar();

            xvelocity = anIntermediateAdiabaticParticle.GetVelocity().GetX();
            yvelocity = anIntermediateAdiabaticParticle.GetVelocity().GetY();
            zvelocity = anIntermediateAdiabaticParticle.GetVelocity().GetZ();

            GammaZ = 1.0/pow(1.0-pow(zvelocity/2.99792e8,2.),0.5);  // fix speed of light.



            LarmorPower = -de/dt*1.602677e-19;

            tLock.unlock();
            fDigitizerCondition.notify_one();  // notify Locust after writing.



             /*
             printf("de is %g and dt is %g and LarmorPower is %g\n", de, dt, LarmorPower);

          	 printf("Kassiopeia says:  tick has happened; continuous time is %g and zvelocity is %f\n", t_poststep, zvelocity);
          	 printf("Kassiopeia says:  fcyc is %g\n", fcyc);
          	 printf("  initial particle momentum: %g\n", anInitialParticle.GetMomentum().Z());
             printf("Mass is %g\n", anInitialParticle.GetMass());
             printf("k.e. = %g eV\n", anInitialParticle.GetKineticEnergy_eV());
             printf("Lorentz factor is %f\n", anInitialParticle.GetLorentzFactor());
             printf("1/(sqrt(1-v^2/c^2) is %f\n", 1.0/pow(1.0-pow(anInitialParticle.GetSpeed()/2.99792e8,2.),0.5));
             printf("FinalParticle().IsActive() is %d\n", aFinalParticle.IsActive());
             */

            delete aDifferentiator;
        	delete anIntegrator;

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

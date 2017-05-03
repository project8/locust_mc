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

    double CyclotronRadiationExtractor::GetGroupVelocity(KSParticle& aFinalParticle)
    {
        double SpeedOfLight = 2.99792458e8; // m/s
        double CutOffFrequency = SpeedOfLight * PI / 10.668e-3; // a in m
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = SpeedOfLight * pow( 1. - pow(CutOffFrequency/(2.*PI*fcyc), 2.) , 0.5);
//        printf("GroupVelocity is %g\n", GroupVelocity); getchar();
    	return GroupVelocity;
    }

    double CyclotronRadiationExtractor::GetCouplingFactor(KSParticle& aFinalParticle)
    {
    	double dim1_wr42 = 10.668e-3; // a in m
    	double vx = aFinalParticle.GetVelocity().GetX();
    	double vy = aFinalParticle.GetVelocity().GetY();
    	double x = aFinalParticle.GetPosition().GetX() + dim1_wr42/2.;
//    	double coupling = fabs(vy)*sin(PI*x/dim1_wr42) / (sin(PI*dim1_wr42/2./dim1_wr42) * pow(vx*vx+vy*vy,0.5));
    	double coupling = 0.63*sin(PI*x/dim1_wr42);  // avg over cyclotron orbit.
    	return coupling*coupling;

    }

    double CyclotronRadiationExtractor::GetDampingFactor(KSParticle& aFinalParticle)
    {
        double fcyc = aFinalParticle.GetCyclotronFrequency();
        double GroupVelocity = GetGroupVelocity(aFinalParticle);
    	double zvelocity = aFinalParticle.GetVelocity().GetZ();
    	double zposition = aFinalParticle.GetPosition().GetZ();
        double GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocity(aFinalParticle),2.),0.5);

    	double fprime_short = fcyc*GammaZ*(1.+zvelocity/GroupVelocity);
    	double phi_short = 2.*PI*2.*(zposition+CENTER_TO_SHORT)/(GroupVelocity/fprime_short);
//        double FieldFromShort = cos(phi_short);  // no resonant enhancement.
        double FieldFromShort = cos(0.) + cos(phi_short); // yes resonant enhancement.
        double CouplingFactor = GetCouplingFactor(aFinalParticle);

        double DampingFactor = CouplingFactor*(1. - FieldFromShort*FieldFromShort);  // can be > 0 or < 0.

//        printf("zposition is %g and lambdaprime is %g\n", zposition, GroupVelocity/fprime_short);
//        printf("CouplingFactor is %f\n", CouplingFactor); getchar();

    	return DampingFactor;
    }



    bool CyclotronRadiationExtractor::ExecutePostStepModifcation( KSParticle& anInitialParticle, KSParticle& aFinalParticle, KSParticleQueue& aQueue )
    {


//    	printf("pre step kinetic energy - 4.84338e-15 is %g\n", anInitialParticle.GetKineticEnergy()- 4.84338e-15); //getchar();
//    	printf("post step kinetic energy - 4.84338e-15 is %g\n", aFinalParticle.GetKineticEnergy()- 4.84338e-15); //getchar();


// adjust power with short.
//    	double DeltaE = GetDampingFactor(aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
//    	aFinalParticle.SetKineticEnergy((aFinalParticle.GetKineticEnergy() - DeltaE));

//    	printf("DeltaE is %g and post fix kinetic energy is %g\n", DeltaE, aFinalParticle.GetKineticEnergy() - 4.84338e-15); getchar();


    	t_poststep = aFinalParticle.GetTime();

        if (t_poststep - t_old > 5.e-10)
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

            GammaZ = 1.0/pow(1.0-pow(zvelocity/GetGroupVelocity(aFinalParticle),2.),0.5);  // speed of light is group velocity

//	    	                fcyc = aFinalParticle.GetCyclotronFrequency();  // inconsistent.  do not use.
					                fcyc = 1.125/dt;
            LarmorPower = -de/dt*1.602677e-19;
            tLock.unlock();
            fDigitizerCondition.notify_one();  // notify Locust after writing.

             t_old = t_poststep;

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

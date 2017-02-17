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

        double X = aFinalParticle.GetPosition().X();
        double Y = aFinalParticle.GetPosition().Y();
        double Z = aFinalParticle.GetPosition().Z();
        double xVelocity = aFinalParticle.GetVelocity().GetX();
        double yVelocity = aFinalParticle.GetVelocity().GetY();
        double zVelocity = aFinalParticle.GetVelocity().GetZ();
        double xMagneticField = aFinalParticle.GetMagneticField().GetX();
        double yMagneticField = aFinalParticle.GetMagneticField().GetY();
        double zMagneticField = aFinalParticle.GetMagneticField().GetZ();
        double mparticle = aFinalParticle.GetMass();
        double qparticle = aFinalParticle.GetCharge();

        locust::ParticleSlim aNewParticle;
        aNewParticle.SetPosition(X,Y,Z);
        aNewParticle.SetVelocityVector(xVelocity,yVelocity,zVelocity);
        aNewParticle.SetMagneticFieldVector(xMagneticField,yMagneticField,zMagneticField);
        aNewParticle.SetMass(mparticle);
        aNewParticle.SetCharge(qparticle);
        aNewParticle.SetTime(t_poststep);
        aNewParticle.SetKinematicProperties();
        
        fNewParticleHistory.push_back(aNewParticle);

        if (t_poststep - t_old > 5.e-10)
        {
        	std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
            tLock.lock();

            de = aFinalParticle.GetKineticEnergy_eV() - anInitialParticle.GetKineticEnergy_eV();
            dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();

            EventModTimeStep = fNewParticleHistory.back().GetTime()-fNewParticleHistory.front().GetTime()+dt;

            //Put in new entries in global ParticleHistory
            fParticleHistory.insert(fParticleHistory.end(),fNewParticleHistory.begin(),fNewParticleHistory.end());
            //Set Spline coefficients
            for(int i=fParticleHistory.size()-fNewParticleHistory.size()-1;i<fParticleHistory.size()-1;i++)
            {
                fParticleHistory[i].SetSpline(fParticleHistory[i+1]);
            }

            fNewParticleHistory.clear();

            //Purge fParticleHistory of overly old entries
            int HistoryMaxSize=5000;
            while(t_poststep-fParticleHistory.front().GetTime()>1e-7 || fParticleHistory.size()>HistoryMaxSize)
            {
                fParticleHistory.pop_front();
            }

            ////Make 100% Sure list is ordered!!!
            //for(int i=0;i<fParticleHistory.size()-1;i++)
            //{
            //    if(fParticleHistory[i].GetTime()>fParticleHistory[i+1].GetTime())
            //    {
            //        printf("Ruh Roh\n");
            //    }
            //}

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

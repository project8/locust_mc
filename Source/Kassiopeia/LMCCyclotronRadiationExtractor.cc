/*
 * LMCCyclotronRadiationExtractor.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCCyclotronRadiationExtractor.hh"
#include "LMCGlobalsDeclaration.hh"
#include "LMCFieldCalculator.hh"

using namespace Kassiopeia;
namespace locust
{

    CyclotronRadiationExtractor::CyclotronRadiationExtractor():
        fP8Phase( 0 ),
        fPitchAngle( -99. )
    {
    }

    CyclotronRadiationExtractor::CyclotronRadiationExtractor( const CyclotronRadiationExtractor& aOrig ):
        fP8Phase( 0 ),
        fPitchAngle( -99. )
    {
    }

    CyclotronRadiationExtractor::~CyclotronRadiationExtractor()
    {
    }

    CyclotronRadiationExtractor* CyclotronRadiationExtractor::Clone() const
    {
        return new CyclotronRadiationExtractor( *this );
    }


    void CyclotronRadiationExtractor::SetP8Phase (int P8Phase )
    {
        fP8Phase = P8Phase;
        Project8Phase = P8Phase;
        if (P8Phase==1)
        {
            CENTER_TO_SHORT = 0.0488; // m, 0.047 is tuned.
            CENTER_TO_ANTENNA = 0.045; // m
        }
        if (P8Phase==2)
        {
            CENTER_TO_SHORT = 0.075; // m
            CENTER_TO_ANTENNA = 0.075; // m  
        }
    }

    bool CyclotronRadiationExtractor::ExecutePreStepModification( KSParticle& anInitialParticle, KSParticleQueue& aQueue )
    {
        return true;
    }


    void CyclotronRadiationExtractor::SetTrajectory( Kassiopeia::KSTrajectory* aTrajectory )
    {
        //this function is being run by KSRoot.cxx at initialization but presently
        //fTrajectory is not staying defined through the stepping.  Suspect binding problem.

        fTrajectory = aTrajectory;

        return;
    }

    locust::Particle CyclotronRadiationExtractor::ExtractKassiopeiaParticle( KSParticle &anInitialParticle, KSParticle &aFinalParticle)
    {
        LMCThreeVector tPosition(aFinalParticle.GetPosition().Components());
        LMCThreeVector tVelocity(aFinalParticle.GetVelocity().Components());
        LMCThreeVector tMagneticField(aFinalParticle.GetMagneticField().Components());
        double tMass = aFinalParticle.GetMass();
        double tCharge = aFinalParticle.GetCharge();
        double tCyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
        double tPitchAngle = aFinalParticle.GetPolarAngleToB();
        double tTime = aFinalParticle.GetTime();


        locust::Particle aNewParticle;
        aNewParticle.SetPosition(tPosition.X(),tPosition.Y(),tPosition.Z());
        aNewParticle.SetVelocityVector(tVelocity.X(),tVelocity.Y(),tVelocity.Z());
        aNewParticle.SetMagneticFieldVector(tMagneticField.X(),tMagneticField.Y(),tMagneticField.Z());
        aNewParticle.SetMass(tMass);
        aNewParticle.SetCharge(tCharge);
        aNewParticle.SetTime(tTime);
        aNewParticle.SetCyclotronFrequency(2.*LMCConst::Pi()*tCyclotronFrequency);
        aNewParticle.SetKinematicProperties();

        if (fPitchAngle == -99.)  // first crossing of center
        {
            if (anInitialParticle.GetPosition().GetZ()/aFinalParticle.GetPosition().GetZ() < 0.)  // trap center
            {
                fPitchAngle = aFinalParticle.GetPolarAngleToB();
                //      	printf("pitch angle is %f\n", fPitchAngle); getchar();

            }
        }
        aNewParticle.SetPitchAngle(fPitchAngle);


        return aNewParticle;

    }


    bool CyclotronRadiationExtractor::ExecutePostStepModification( KSParticle& anInitialParticle, KSParticle& aFinalParticle, KSParticleQueue& aQueue )
    {

        FieldCalculator aFieldCalculator;
        double DeltaE=0.;

        //		       printf("fcyc is %g\n", aFinalParticle.GetCyclotronFrequency()); getchar();

        //	printf("dE/dt is %g\n", (aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy())/(aFinalParticle.GetTime() - anInitialParticle.GetTime())); getchar();

        if(fP8Phase==1)
        {
            DeltaE = aFieldCalculator.GetDampingFactorPhase1(anInitialParticle, aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            aFinalParticle.SetKineticEnergy((anInitialParticle.GetKineticEnergy() + DeltaE));
        }
        if(fP8Phase==2)
        {
            DeltaE = aFieldCalculator.GetDampingFactorPhase2(anInitialParticle, aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            aFinalParticle.SetKineticEnergy((anInitialParticle.GetKineticEnergy() + DeltaE));
        }

        if (!fDoneWithSignalGeneration)  // if Locust is still acquiring voltages.
        {

            if (t_old == 0.) 
            {
                fPitchAngle = -99.;  // new electron needs central pitch angle reset.
                if (fParticleHistory.size()) fParticleHistory.clear();
            }
            double t_poststep = aFinalParticle.GetTime();
            fNewParticleHistory.push_back(ExtractKassiopeiaParticle(anInitialParticle, aFinalParticle));

            if (t_poststep - t_old >= fDigitizerTimeStep) //take a digitizer sample every 5e-10s
            {
                std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
                tLock.lock();

                int tHistoryMaxSize;

                //Dont want to check .back() of history if it is empty! -> Segfault
                if(fParticleHistory.size() && (fNewParticleHistory.back().GetTime() < fParticleHistory.back().GetTime()))
                {
                    t_poststep = 0.;
                    fParticleHistory.clear();
                }


                //Put in new entries in global ParticleHistory
                fParticleHistory.insert(fParticleHistory.end(),fNewParticleHistory.begin(),fNewParticleHistory.end());

                for(int i=fParticleHistory.size()-fNewParticleHistory.size()-1;i<fParticleHistory.size()-1;i++)
                {
                    fParticleHistory[i].SetSpline(fParticleHistory[i+1]);
                }

                tHistoryMaxSize = 5000;

                fNewParticleHistory.clear();

                //Purge fParticleHistory of overly old entries	    
                while(t_poststep-fParticleHistory.front().GetTime()>1e-7 || fParticleHistory.size() > tHistoryMaxSize)
                {
                    fParticleHistory.pop_front();
                }

                tLock.unlock();
                fDigitizerCondition.notify_one();  // notify Locust after writing.

            }
        } // DoneWithSignalGeneration

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

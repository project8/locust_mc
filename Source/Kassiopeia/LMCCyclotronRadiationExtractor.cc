#include "LMCCyclotronRadiationExtractor.hh"
#include "LMCGlobalsDeclaration.hh"
#include "KSModifiersMessage.h"
#include "LMCFieldCalculator.hh"

namespace locust
{
    CyclotronRadiationExtractor::CyclotronRadiationExtractor() :
            fP8Phase( 0 ),
            fNewParticleHistory(),
            fPitchAngle( -99. )
    {
    }

    CyclotronRadiationExtractor::CyclotronRadiationExtractor(const CyclotronRadiationExtractor &aCopy) : KSComponent(),
            fP8Phase( aCopy.fP8Phase ),
            fNewParticleHistory(),
            fPitchAngle( aCopy.fPitchAngle )
    {
    }
    CyclotronRadiationExtractor* CyclotronRadiationExtractor::Clone() const
    {
        return new CyclotronRadiationExtractor( *this );
    }
    CyclotronRadiationExtractor::~CyclotronRadiationExtractor()
    {
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


    locust::Particle CyclotronRadiationExtractor::ExtractKassiopeiaParticle( Kassiopeia::KSParticle &anInitialParticle, Kassiopeia::KSParticle &aFinalParticle)
    {
        LMCThreeVector tPosition(aFinalParticle.GetPosition().Components());
        LMCThreeVector tVelocity(aFinalParticle.GetVelocity().Components());
        LMCThreeVector tMagneticField(aFinalParticle.GetMagneticField().Components());
        double tMass = aFinalParticle.GetMass();
        double tCharge = aFinalParticle.GetCharge();
        double tCyclotronFrequency = aFinalParticle.GetCyclotronFrequency();
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
            }
        }
        aNewParticle.SetPitchAngle(fPitchAngle);

        return aNewParticle;

    }



    bool CyclotronRadiationExtractor::ExecutePostStepModification( Kassiopeia::KSParticle& anInitialParticle,
                                                                   Kassiopeia::KSParticle& aFinalParticle,
                                                                   Kassiopeia::KSParticleQueue& aParticleQueue )
    {

        FieldCalculator aFieldCalculator;
        double DeltaE=0.;

        //		       printf("fcyc is %g\n", aFinalParticle.GetCyclotronFrequency()); getchar();

        //	printf("dE/dt is %g\n", (aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy())/(aFinalParticle.GetTime() - anInitialParticle.GetTime())); getchar();

        if(fP8Phase==1)
        {
            DeltaE = aFieldCalculator.GetDampingFactorPhase1(aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            aFinalParticle.SetKineticEnergy((anInitialParticle.GetKineticEnergy() + DeltaE));
        }
        if(fP8Phase==2)  // this code can be commented out to save time as DeltaE will be small.
        {
            DeltaE = aFieldCalculator.GetDampingFactorPhase2(aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            aFinalParticle.SetKineticEnergy((anInitialParticle.GetKineticEnergy() + DeltaE));
        }

        if (!fDoneWithSignalGeneration)  // if Locust is still acquiring voltages.
        {

            if (t_old == 0.) 
            {
                fPitchAngle = -99.;  // new electron needs central pitch angle reset.
            }
            double t_poststep = aFinalParticle.GetTime();
            fNewParticleHistory.push_back(ExtractKassiopeiaParticle(anInitialParticle, aFinalParticle));

            if (t_poststep - t_old >= fKassTimeStep) //take a digitizer sample every KassTimeStep
            {
                std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
                tLock.lock();

                unsigned tHistoryMaxSize;

                //Dont want to check .back() of history if it is empty! -> Segfault
                if(fParticleHistory.size() && (fNewParticleHistory.back().GetTime() < fParticleHistory.back().GetTime()))
                {
                    t_poststep = 0.;
                    fParticleHistory.clear();
                }
                //                else {printf("history is empty at t_old %g\n", t_old); getchar();}


                //Put in new entries in global ParticleHistory
                fParticleHistory.insert(fParticleHistory.end(),fNewParticleHistory.begin(),fNewParticleHistory.end());

                for(unsigned i=fParticleHistory.size()-fNewParticleHistory.size()-1;i<fParticleHistory.size()-1;i++)
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

        return false;
    }

/*
    STATICINT sKSRootModifierDict =
            Kassiopeia::KSDictionary< CyclotronRadiationExtractor >::AddCommand( &CyclotronRadiationExtractor::AddModifier,
                    &CyclotronRadiationExtractor::RemoveModifier,
                    "add_modifier",
                    "remove_modifier" );
*/
}

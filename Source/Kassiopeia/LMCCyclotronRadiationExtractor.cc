#include "LMCCyclotronRadiationExtractor.hh"
#include "KSModifiersMessage.h"


namespace locust
{

    LOGGER( lmclog, "CyclotronRadiationExtractor" );

    CyclotronRadiationExtractor::CyclotronRadiationExtractor() :
            fNewParticleHistory(),
			fFieldCalculator( NULL ),
            fPitchAngle( -99. ),
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    	Configure();
    }

    CyclotronRadiationExtractor::CyclotronRadiationExtractor(const CyclotronRadiationExtractor &aCopy) : KSComponent(),
            fNewParticleHistory(),
			fFieldCalculator( NULL ),
            fPitchAngle( aCopy.fPitchAngle ),
            fInterface( aCopy.fInterface )
    {
    	Configure();
    }

    CyclotronRadiationExtractor* CyclotronRadiationExtractor::Clone() const
    {
        return new CyclotronRadiationExtractor( *this );
    }
    CyclotronRadiationExtractor::~CyclotronRadiationExtractor()
    {
    }

    bool CyclotronRadiationExtractor::Configure()
    {
    	fFieldCalculator = new FieldCalculator();
    	if(!fFieldCalculator->ConfigureByInterface())
    	{
    	   LERROR(lmclog,"Error configuring receiver FieldCalculator class from CyclotronRadiationExtractor.");
    	   exit(-1);
    	}
    	return true;
    }


    void CyclotronRadiationExtractor::SetP8Phase (int P8Phase )
    {
        fInterface->fProject8Phase = P8Phase;
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

    bool CyclotronRadiationExtractor::ExecutePreStepModification( Kassiopeia::KSParticle& , Kassiopeia::KSParticleQueue&  )
    {
        return true;
    }

    bool CyclotronRadiationExtractor::ExecutePostStepModification( Kassiopeia::KSParticle& anInitialParticle,
                                                                   Kassiopeia::KSParticle& aFinalParticle,
                                                                   Kassiopeia::KSParticleQueue& aParticleQueue )
    {

        double DeltaE=0.;

        if(fInterface->fProject8Phase==1)
        {
            DeltaE = fFieldCalculator->GetDampingFactorPhase1(aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            aFinalParticle.SetKineticEnergy((anInitialParticle.GetKineticEnergy() + DeltaE));
        }
        if(fInterface->fProject8Phase==2)
        {
            DeltaE = fFieldCalculator->GetDampingFactorPhase2(aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            aFinalParticle.SetKineticEnergy((anInitialParticle.GetKineticEnergy() + DeltaE));
        }

        if(fInterface->fProject8Phase==4) // Cavity
        {
            if (fInterface->fE_Gun)
            {
            	DeltaE = fFieldCalculator->GetDampingFactorPhase1(aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            }
            else
            {
//		std::cout << "About to call GetDampingFactor Cavity from LMCCyclotronRadiationExtractor" << std::endl;
            	DeltaE = fFieldCalculator->GetDampingFactorCavity(1, 1, 1, aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy()); //Infrastructure needs to be generalized to other modes
//		std::cout << "DeltaE is " << DeltaE << std::endl;
            }
            if (fInterface->fBackReaction)
            {
            	aFinalParticle.SetKineticEnergy((anInitialParticle.GetKineticEnergy() + DeltaE));
            }
            else
            {
            	// Do not apply any power correction, even though it was calculated above.
            }
        }



        if (!fInterface->fDoneWithSignalGeneration)  // if Locust is still acquiring voltages.
        {
            if (fInterface->fTOld == 0.)
            {
            	fPitchAngle = -99.;  // new electron needs central pitch angle reset.
            	double dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();
//		std::cout << "Setting NFilterBinsRequired in CyclotronRadiationExtractor" << std::endl;
                fFieldCalculator->SetNFilterBinsRequired( dt );
            }

            double t_poststep = aFinalParticle.GetTime();
            fNewParticleHistory.push_back(ExtractKassiopeiaParticle(anInitialParticle, aFinalParticle));

            if (t_poststep - fInterface->fTOld >= fInterface->fKassTimeStep) //take a digitizer sample every KassTimeStep
            {

                std::unique_lock< std::mutex >tLock( fInterface->fMutexDigitizer, std::defer_lock );  // lock access to mutex before writing to globals.
                tLock.lock();


                unsigned tHistoryMaxSize;

                //Dont want to check .back() of history if it is empty! -> Segfault
                if(fInterface->fParticleHistory.size() && (fNewParticleHistory.back().GetTime() < fInterface->fParticleHistory.back().GetTime()))
                {
                    t_poststep = 0.;
                    fInterface->fParticleHistory.clear();
                }
                //                else {printf("history is empty at t_old %g\n", t_old); getchar();}


                //Put in new entries in global ParticleHistory
                fInterface->fParticleHistory.insert(fInterface->fParticleHistory.end(),fNewParticleHistory.begin(),fNewParticleHistory.end());

                for(unsigned i=fInterface->fParticleHistory.size()-fNewParticleHistory.size()-1;i<fInterface->fParticleHistory.size()-1;i++)
                {
                    fInterface->fParticleHistory[i].SetSpline(fInterface->fParticleHistory[i+1]);
                }

                tHistoryMaxSize = 5000;

                fNewParticleHistory.clear();

                //Purge fParticleHistory of overly old entries	    
                while(t_poststep-fInterface->fParticleHistory.front().GetTime()>1e-7 || fInterface->fParticleHistory.size() > tHistoryMaxSize)
                {
                    fInterface->fParticleHistory.pop_front();
                }

                tLock.unlock();

                fInterface->fDigitizerCondition.notify_one();  // notify Locust after writing.

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

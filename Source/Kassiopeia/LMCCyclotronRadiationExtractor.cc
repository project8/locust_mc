#include "LMCCyclotronRadiationExtractor.hh"
#include "LMCGlobalsDeclaration.hh"
#include "KSModifiersMessage.h"
#include "LMCFieldCalculator.hh"

namespace locust
{
    CyclotronRadiationExtractor::CyclotronRadiationExtractor() :
        fModifiers(128),
        fModifier( NULL ),
        fInitialParticle( NULL ),
        fFinalParticle( NULL ),
        fParticleQueue( NULL ),
        fP8Phase( 0 ),
        fPitchAngle( -99. )
    {
    }

    CyclotronRadiationExtractor::CyclotronRadiationExtractor(const CyclotronRadiationExtractor &aCopy) : KSComponent(),
        fModifiers( aCopy.fModifiers ),
        fModifier( aCopy.fModifier),
        fInitialParticle( aCopy.fInitialParticle ),
        fFinalParticle( aCopy.fFinalParticle ),
        fParticleQueue( aCopy.fParticleQueue )
    {
    }
    CyclotronRadiationExtractor* CyclotronRadiationExtractor::Clone() const
    {
        return new CyclotronRadiationExtractor( *this );
    }
    CyclotronRadiationExtractor::~CyclotronRadiationExtractor()
    {
    }

    void CyclotronRadiationExtractor::AddModifier(KSStepModifier *aModifier)
    {
        fModifiers.AddElement( aModifier );
        return;
    }
    void CyclotronRadiationExtractor::RemoveModifier(KSStepModifier *aModifier)
    {
        fModifiers.RemoveElement( aModifier );
        return;
    }
  void CyclotronRadiationExtractor::SetStep( Kassiopeia::KSStep* aStep )
    {
        fStep = aStep;
        fInitialParticle = &(aStep->InitialParticle());
        fModifierParticle = &(aStep->TerminatorParticle());
        fFinalParticle = &(aStep->FinalParticle());
        fParticleQueue = &(aStep->ParticleQueue());

        return;
    }

    void CyclotronRadiationExtractor::PushUpdateComponent()
    {
        for( int tIndex = 0; tIndex < fModifiers.End(); tIndex++ )
        {
            fModifiers.ElementAt( tIndex )->PushUpdate();
        }
    }

    void CyclotronRadiationExtractor::PushDeupdateComponent()
    {
        for( int tIndex = 0; tIndex < fModifiers.End(); tIndex++ )
        {
            fModifiers.ElementAt( tIndex )->PushDeupdate();
        }
    }

    bool CyclotronRadiationExtractor::ExecutePreStepModification()
    {
        //the following disables any change made to the initial particle, why?
        *fModifierParticle = *fInitialParticle;
        fStep->ModifierName().clear();
        fStep->ModifierFlag() = false;

        if( fModifiers.End() == 0 )
        {
            modmsg_debug( "modifier calculation:" << eom )
            modmsg_debug( "  no modifier active" << eom )
            modmsg_debug( "  modifier name: <" << fStep->GetModifierName() << ">" << eom )
            modmsg_debug( "  modifier flag: <" << fStep->GetModifierFlag() << ">" << eom )

            modmsg_debug( "modifier calculation terminator particle state: " << eom )
            modmsg_debug( "  modifier particle space: <" << (fModifierParticle->GetCurrentSpace() ? fModifierParticle->GetCurrentSpace()->GetName() : "" ) << ">" << eom )
            modmsg_debug( "  modifier particle surface: <" << (fModifierParticle->GetCurrentSurface() ? fModifierParticle->GetCurrentSurface()->GetName() : "" ) << ">" << eom )
            modmsg_debug( "  modifier particle time: <" << fModifierParticle->GetTime() << ">" << eom )
            modmsg_debug( "  modifier particle length: <" << fModifierParticle->GetLength() << ">" << eom )
            modmsg_debug( "  modifier particle position: <" << fModifierParticle->GetPosition().X() << ", " << fModifierParticle->GetPosition().Y() << ", " << fModifierParticle->GetPosition().Z() << ">" << eom )
            modmsg_debug( "  modifier particle momentum: <" << fModifierParticle->GetMomentum().X() << ", " << fModifierParticle->GetMomentum().Y() << ", " << fModifierParticle->GetMomentum().Z() << ">" << eom )
            modmsg_debug( "  modifier particle kinetic energy: <" << fModifierParticle->GetKineticEnergy_eV() << ">" << eom )
            modmsg_debug( "  modifier particle electric field: <" << fModifierParticle->GetElectricField().X() << "," << fModifierParticle->GetElectricField().Y() << "," << fModifierParticle->GetElectricField().Z() << ">" << eom )
            modmsg_debug( "  modifier particle magnetic field: <" << fModifierParticle->GetMagneticField().X() << "," << fModifierParticle->GetMagneticField().Y() << "," << fModifierParticle->GetMagneticField().Z() << ">" << eom )
            modmsg_debug( "  modifier particle angle to magnetic field: <" << fModifierParticle->GetPolarAngleToB() << ">" << eom )

            return false; //changes to inital particle state disabled
        }

        fStep->ModifierFlag() = ExecutePreStepModification( *fModifierParticle, *fParticleQueue );

        (void) fStep->ModifierFlag();
        //hasChangedState is unused because we are operating on the modifier particle
        //this disables any changes to the intial particle

        if( fStep->ModifierFlag() == true )
        {
            modmsg_debug( "modifier calculation:" << eom )
            modmsg_debug( "  modification may occur" << eom )
        }
        else
        {
            modmsg_debug( "modifier calculation:" << eom )
            modmsg_debug( "  modification will not occur" << eom )
        }

        modmsg_debug( "modifier calculation modifier particle state: " << eom )
        modmsg_debug( "  modifier particle space: <" << (fModifierParticle->GetCurrentSpace() ? fModifierParticle->GetCurrentSpace()->GetName() : "" ) << ">" << eom )
        modmsg_debug( "  modifier particle surface: <" << (fModifierParticle->GetCurrentSurface() ? fModifierParticle->GetCurrentSurface()->GetName() : "" ) << ">" << eom )
        modmsg_debug( "  modifier particle time: <" << fModifierParticle->GetTime() << ">" << eom )
        modmsg_debug( "  modifier particle length: <" << fModifierParticle->GetLength() << ">" << eom )
        modmsg_debug( "  modifier particle position: <" << fModifierParticle->GetPosition().X() << ", " << fModifierParticle->GetPosition().Y() << ", " << fModifierParticle->GetPosition().Z() << ">" << eom )
        modmsg_debug( "  modifier particle momentum: <" << fModifierParticle->GetMomentum().X() << ", " << fModifierParticle->GetMomentum().Y() << ", " << fModifierParticle->GetMomentum().Z() << ">" << eom )
        modmsg_debug( "  modifier particle kinetic energy: <" << fModifierParticle->GetKineticEnergy_eV() << ">" << eom )
        modmsg_debug( "  modifier particle electric field: <" << fModifierParticle->GetElectricField().X() << "," << fModifierParticle->GetElectricField().Y() << "," << fModifierParticle->GetElectricField().Z() << ">" << eom )
        modmsg_debug( "  modifier particle magnetic field: <" << fModifierParticle->GetMagneticField().X() << "," << fModifierParticle->GetMagneticField().Y() << "," << fModifierParticle->GetMagneticField().Z() << ">" << eom )
        modmsg_debug( "  modifier particle angle to magnetic field: <" << fModifierParticle->GetPolarAngleToB() << ">" << eom )

        return false; //changes to initial particle state disabled
    }

    void CyclotronRadiationExtractor::SetTrajectory( Kassiopeia::KSTrajectory* aTrajectory )
    {
        //this function is being run by KSRoot.cxx at initialization but presently
        //fTrajectory is not staying defined through the stepping.  Suspect binding problem.

        fTrajectory = aTrajectory;

        return;
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
                //      	printf("pitch angle is %f\n", fPitchAngle); getchar();
            }
        }
        aNewParticle.SetPitchAngle(fPitchAngle);

        return aNewParticle;

    }



    bool CyclotronRadiationExtractor::ExecutePostStepModification()
    {


        bool hasChangedState = ExecutePostStepModification( *fModifierParticle, *fFinalParticle, *fParticleQueue );
        fFinalParticle->ReleaseLabel( fStep->ModifierName() );

        modmsg_debug( "modifier execution:" << eom )
        modmsg_debug( "  terminator name: <" << fStep->TerminatorName() << ">" << eom )
        modmsg_debug( "  step continuous time: <" << fStep->ContinuousTime() << ">" << eom )
        modmsg_debug( "  step continuous length: <" << fStep->ContinuousLength() << ">" << eom )
        modmsg_debug( "  step continuous energy change: <" << fStep->ContinuousEnergyChange() << ">" << eom )
        modmsg_debug( "  step continuous momentum change: <" << fStep->ContinuousMomentumChange() << ">" << eom )
        modmsg_debug( "  step discrete secondaries: <" << fStep->DiscreteSecondaries() << ">" << eom )
        modmsg_debug( "  step discrete energy change: <" << fStep->DiscreteEnergyChange() << ">" << eom )
        modmsg_debug( "  step discrete momentum change: <" << fStep->DiscreteMomentumChange() << ">" << eom )

        modmsg_debug( "modifier final particle state: " << eom )
        modmsg_debug( "  final particle space: <" << (fModifierParticle->GetCurrentSpace() ? fModifierParticle->GetCurrentSpace()->GetName() : "" ) << ">" << eom )
        modmsg_debug( "  final particle surface: <" << (fModifierParticle->GetCurrentSurface() ? fModifierParticle->GetCurrentSurface()->GetName() : "" ) << ">" << eom )
        modmsg_debug( "  final particle time: <" << fModifierParticle->GetTime() << ">" << eom )
        modmsg_debug( "  final particle length: <" << fModifierParticle->GetLength() << ">" << eom )
        modmsg_debug( "  final particle position: <" << fModifierParticle->GetPosition().X() << ", " << fModifierParticle->GetPosition().Y() << ", " << fModifierParticle->GetPosition().Z() << ">" << eom )
        modmsg_debug( "  final particle momentum: <" << fModifierParticle->GetMomentum().X() << ", " << fModifierParticle->GetMomentum().Y() << ", " << fModifierParticle->GetMomentum().Z() << ">" << eom )
        modmsg_debug( "  final particle kinetic energy: <" << fModifierParticle->GetKineticEnergy_eV() << ">" << eom )
        modmsg_debug( "  final particle electric field: <" << fModifierParticle->GetElectricField().X() << "," << fModifierParticle->GetElectricField().Y() << "," << fModifierParticle->GetElectricField().Z() << ">" << eom )
        modmsg_debug( "  final particle magnetic field: <" << fModifierParticle->GetMagneticField().X() << "," << fModifierParticle->GetMagneticField().Y() << "," << fModifierParticle->GetMagneticField().Z() << ">" << eom )
        modmsg_debug( "  final particle angle to magnetic field: <" << fModifierParticle->GetPolarAngleToB() << ">" << eom )

        return hasChangedState;
    }

  bool CyclotronRadiationExtractor::ExecutePreStepModification( Kassiopeia::KSParticle& anInitialParticle,
								Kassiopeia::KSParticleQueue& aParticleQueue )
    {
        bool hasChangedState = false;
        for( int tIndex = 0; tIndex < fModifiers.End(); tIndex++ )
        {
            bool changed = fModifiers.ElementAt( tIndex )->ExecutePreStepModification( anInitialParticle, aParticleQueue );
            if(changed){hasChangedState = true;};
        }

        return hasChangedState;
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
        if(fP8Phase==2)  // this code be commented out to save time as DeltaE will be small.
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

            if (t_poststep - t_old >= fDigitizerTimeStep) //take a digitizer sample every 5e-10s
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


        bool hasChangedState = false;
        for( int tIndex = 0; tIndex < fModifiers.End(); tIndex++ )
        {
            bool changed = fModifiers.ElementAt( tIndex )->ExecutePostStepModification( anInitialParticle,
                                                                        aFinalParticle,
                                                                        aParticleQueue );
            if(changed){hasChangedState = true;};
        }

        return hasChangedState;
    }


    STATICINT sKSRootModifierDict =
														Kassiopeia::KSDictionary< CyclotronRadiationExtractor >::AddCommand( &CyclotronRadiationExtractor::AddModifier,
                                                            &CyclotronRadiationExtractor::RemoveModifier,
                                                            "add_modifier",
                                                            "remove_modifier" );

}

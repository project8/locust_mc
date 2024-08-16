#include "LMCCyclotronRadiationExtractor.hh"
#include "KSModifiersMessage.h"
#include <chrono>
#include <thread>


namespace locust
{

    LOGGER( lmclog, "CyclotronRadiationExtractor" );

    CyclotronRadiationExtractor::CyclotronRadiationExtractor() :
            fNewParticleHistory(),
            fFieldCalculator( NULL ),
            fPitchAngle( -99. ),
            fLMCTrackID( -2 ),
            fT0trapMin( 0. ),
            fNCrossings( 0 ),
            fSampleIndex( 0 ),
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    CyclotronRadiationExtractor::CyclotronRadiationExtractor(const CyclotronRadiationExtractor &aCopy) : KSComponent(),
            fNewParticleHistory(),
            fFieldCalculator( NULL ),
            fPitchAngle( aCopy.fPitchAngle ),
            fLMCTrackID( aCopy.fLMCTrackID ),
            fT0trapMin( aCopy.fT0trapMin ),
            fNCrossings( aCopy.fNCrossings ),
            fSampleIndex( aCopy.fSampleIndex ),
            fInterface( aCopy.fInterface )
    {
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
        if (fInterface->fProject8Phase > 0)
        {
            if(!fFieldCalculator->ConfigureByInterface())
            {
                LERROR(lmclog,"Error configuring receiver FieldCalculator class from CyclotronRadiationExtractor.");
                exit(-1);
            }
        }
        return true;
    }


    void CyclotronRadiationExtractor::SetP8Phase (int P8Phase )
    {
        fInterface->fProject8Phase = P8Phase;
        Configure();
    }

    bool CyclotronRadiationExtractor::UpdateTrackProperties( Kassiopeia::KSParticle &aFinalParticle, unsigned index, bool bStart )
    {
    	double tTime = index / fInterface->aRunParameter->fSamplingRateMHz / 1.e6 / fInterface->aRunParameter->fDecimationFactor;
#ifdef ROOT_FOUND
    	if (bStart)
    	{
            fStartingIndex = index;
            fInterface->aTrack.StartTime = tTime;
            double tX = aFinalParticle.GetPosition().X();
            double tY = aFinalParticle.GetPosition().Y();
            fInterface->aTrack.Radius = pow(tX*tX + tY*tY, 0.5);
            fInterface->aTrack.RadialPhase = calcOrbitPhase(tX, tY);
            fInterface->aTrack.StartingEnergy_eV = LMCConst::kB_eV() / LMCConst::kB() * aFinalParticle.GetKineticEnergy();
    	}
    	else
    	{
            fInterface->aTrack.EndTime = tTime;
            unsigned nElapsedSamples = index - fStartingIndex;
            fInterface->aTrack.AvgFrequency = ( fInterface->aTrack.AvgFrequency * nElapsedSamples + aFinalParticle.GetCyclotronFrequency() ) / ( nElapsedSamples + 1);
            fInterface->aTrack.TrackLength = tTime - fInterface->aTrack.StartTime;
            fInterface->aTrack.Slope = (fInterface->aTrack.EndFrequency - fInterface->aTrack.StartFrequency) / (fInterface->aTrack.TrackLength);
    	}
#endif

        return true;
    }

    double CyclotronRadiationExtractor::calcOrbitPhase(double tX, double tY)
    {
    	double phase = 0.;
        if ((fabs(tX) > 0.))
    	{
    		phase = atan(tY/tX);
    	}

    	phase += quadrantOrbitCorrection(phase, tY);
    	return phase;
    }

    double CyclotronRadiationExtractor::quadrantOrbitCorrection(double phase, double tY)
    {
    	double phaseCorrection = 0.;
    	if (((phase < 0.)&&(tY > 0.)) || ((phase > 0.)&&(tY < 0.)))
    		phaseCorrection = LMCConst::Pi();

    	return phaseCorrection;
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


        if (anInitialParticle.GetPosition().GetZ()/aFinalParticle.GetPosition().GetZ() < 0.)  // trap center
        {
            fNCrossings += 1;
            if (fPitchAngle == -99.)  // first crossing of center
            {
                fPitchAngle = aFinalParticle.GetPolarAngleToB();
                fT0trapMin = aFinalParticle.GetTime();
#ifdef ROOT_FOUND
                fInterface->aTrack.PitchAngle = aFinalParticle.GetPolarAngleToB();
                fInterface->aTrack.StartFrequency = aFinalParticle.GetCyclotronFrequency();
                double tLOfrequency = fInterface->aRunParameter->fLOfrequency; // Hz
                double tSamplingRate = fInterface->aRunParameter->fSamplingRateMHz; // MHz
                fInterface->aTrack.LOFrequency = tLOfrequency;
                fInterface->aTrack.RandomSeed = fInterface->aRunParameter->fRandomSeed;
                fInterface->aTrack.OutputStartFrequency = fInterface->aTrack.StartFrequency - tLOfrequency + tSamplingRate * 1.e6 / 2.;
#endif
            }
            else
            {
#ifdef ROOT_FOUND
                fInterface->aTrack.EndFrequency = aFinalParticle.GetCyclotronFrequency();
                fInterface->aTrack.AvgAxialFrequency = fNCrossings / 2. / ( aFinalParticle.GetTime() - fT0trapMin );
#endif
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
            if (fInterface->fbWaveguide)
            {
            	DeltaE = fFieldCalculator->GetDampingFactorPhase1(aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
            }
            else
            {
            	DeltaE = fFieldCalculator->GetDampingFactorCavity(aFinalParticle)*(aFinalParticle.GetKineticEnergy() - anInitialParticle.GetKineticEnergy());
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

            if ( aFinalParticle.GetIndexNumber() > fLMCTrackID ) // check for new track
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1000));
                fLMCTrackID = aFinalParticle.GetIndexNumber();
                fPitchAngle = -99.;  // new electron needs central pitch angle reset.
                double dt = aFinalParticle.GetTime() - anInitialParticle.GetTime();
                fFieldCalculator->SetNFilterBinsRequired( dt );
                UpdateTrackProperties( aFinalParticle, fInterface->fSampleIndex, 1 );
                LPROG(lmclog,"Updated recorded track properties at sample " << fInterface->fSampleIndex );
            }


            double t_poststep = aFinalParticle.GetTime();
            fNewParticleHistory.push_back(ExtractKassiopeiaParticle(anInitialParticle, aFinalParticle));

            if (t_poststep - fInterface->fTOld >= fInterface->fKassTimeStep) //take a digitizer sample every KassTimeStep
            {

                fSampleIndex = fInterface->fSampleIndex; // record Locust sample index before locking
                UpdateTrackProperties( aFinalParticle, fSampleIndex, 0 );  // Keep recording the track candidate end time.

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

                int tTriggerConfirm = 0;

                while ( !(fSampleIndex < fInterface->fSampleIndex) && (tTriggerConfirm < fInterface->fTriggerConfirm) )
                {
                    // If the Locust sample index has not advanced yet, keep checking it.
                    tTriggerConfirm += 1;
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    if ( tTriggerConfirm % 1000 == 0 )
                    {
                    	LPROG(lmclog,"Checking the digitizer synchronization, tTriggerConfirm index = " << tTriggerConfirm );
                    }

                    if ( ( tTriggerConfirm > fInterface->fTriggerConfirm - 3) && ( fSampleIndex < fInterface->fFastRecordLength-1 ) )
                    {
                        LPROG(lmclog,"Checking the digitizer synchronization, tTriggerConfirm index = " << tTriggerConfirm);
                        LPROG(lmclog,"Checking the digitizer synchronization, at fast sample = " << fSampleIndex);
                        LPROG(lmclog,"Checking the digitizer synchronization, at Locust fast sample = " << fInterface->fSampleIndex);
                        LPROG(lmclog,"Fast record length = " << fInterface->fFastRecordLength);
                        std::this_thread::sleep_for(std::chrono::milliseconds(10000));
                        if ( !(fSampleIndex < fInterface->fSampleIndex) )
                        {
                            LPROG(lmclog,"Checking the digitizer synchronization again.  ");
                            std::this_thread::sleep_for(std::chrono::milliseconds(10000));
                            if ( !(fSampleIndex < fInterface->fSampleIndex) )
                            {
                                LERROR(lmclog,"Locust digitizer sample index has not advanced properly.  "
                         	    		"Please either resubmit the job, or check HPC status.");
                                LERROR(lmclog, "tTriggerConfirm, fSampleIndex are " << tTriggerConfirm << " and " << fSampleIndex);
                                exit(-1);  // TO-DO:  throw this exception to be caught properly by scarab, as in LocustSim.cc .
                            }
                        }
                    }
                }

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

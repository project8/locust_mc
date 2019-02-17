
#ifndef LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_
#define LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_

#include "KSStepModifier.h"
#include "KSStep.h"
#include "KSList.h"
#include "KSParticle.h"
#include "KSTrajectory.h"
#include "LMCParticle.hh"
#include "KToolbox.h"



namespace locust
{

    class KSTrack;

    class CyclotronRadiationExtractor :
    public Kassiopeia::KSComponentTemplate< CyclotronRadiationExtractor, Kassiopeia::KSStepModifier >
    {
    public:
        CyclotronRadiationExtractor();
        CyclotronRadiationExtractor( const CyclotronRadiationExtractor& aCopy );
        CyclotronRadiationExtractor* Clone() const;
        virtual ~CyclotronRadiationExtractor();

        //**********
        // modifier
        //**********

    public:
      bool ExecutePreStepModification( Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticleQueue& aQueue );
        bool ExecutePostStepModification( Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle, Kassiopeia::KSParticleQueue& aQueue );

        //***********
        //composition
        //***********

    public:
        void AddModifier( Kassiopeia::KSStepModifier* aModifier );
        void RemoveModifier( Kassiopeia::KSStepModifier* aModifier );

    private:
      Kassiopeia::KSList< Kassiopeia::KSStepModifier > fModifiers;
      Kassiopeia::KSStepModifier* fModifier;

        //******
        //action
        //******

    public:
      void SetStep( Kassiopeia::KSStep* aStep );

        bool ExecutePreStepModification();
        bool ExecutePostStepModification();


        locust::Particle ExtractKassiopeiaParticle( Kassiopeia::KSParticle &anInitialParticle, Kassiopeia::KSParticle &aFinalParticle);

            void SetTrajectory( Kassiopeia::KSTrajectory* aTrajectory );
            void SetP8Phase( int P8Phase );


        virtual void PushUpdateComponent();
        virtual void PushDeupdateComponent();

    private:
      Kassiopeia::KSStep* fStep;
        const Kassiopeia::KSParticle* fInitialParticle;
        Kassiopeia::KSParticle* fModifierParticle;
        Kassiopeia::KSParticle* fFinalParticle;
        Kassiopeia::KSParticleQueue* fParticleQueue;
        Kassiopeia::KSTrajectory* fTrajectory;
        int fP8Phase; // 1, 2, 3, or 4.
        std::deque<locust::Particle> fNewParticleHistory;
        double fPitchAngle;
        double fCentralPower;
        double fCentralZVelocity;

    };


}

#endif

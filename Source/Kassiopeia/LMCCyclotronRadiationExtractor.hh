
#ifndef LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_
#define LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_

#include "KSStepModifier.h"
#include "KSTrajectory.h"

#include "LMCParticle.hh"

#include <deque>

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

            locust::Particle ExtractKassiopeiaParticle( Kassiopeia::KSParticle &anInitialParticle, Kassiopeia::KSParticle &aFinalParticle);

            void SetTrajectory( Kassiopeia::KSTrajectory* aTrajectory );
            void SetP8Phase( int P8Phase );


            virtual void PushUpdateComponent();
            virtual void PushDeupdateComponent();

        private:
            int fP8Phase; // 1, 2, 3, or 4.
            std::deque<locust::Particle> fNewParticleHistory;
            double fPitchAngle;

    };


}

#endif

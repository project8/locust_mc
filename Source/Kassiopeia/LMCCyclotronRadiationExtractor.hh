
#ifndef LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_
#define LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_

#include "KSStepModifier.h"
#include "KSTrajectory.h"

#include "LMCFieldCalculator.hh"
#include "LMCKassLocustInterface.hh"
#include "LMCParticle.hh"
#include "LMCException.hh"

#ifdef ROOT_FOUND
    #include "LMCRootTreeWriter.hh"
#endif



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
            bool Configure();


            //**********
            // modifier
            //**********

        public:
            bool ExecutePreStepModification( Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticleQueue& aQueue );
            bool ExecutePostStepModification( Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle, Kassiopeia::KSParticleQueue& aQueue );

            locust::Particle ExtractKassiopeiaParticle( Kassiopeia::KSParticle &anInitialParticle, Kassiopeia::KSParticle &aFinalParticle);


            void SetTrajectory( Kassiopeia::KSTrajectory* aTrajectory );
            void SetP8Phase( int P8Phase );
            bool UpdateTrackProperties( Kassiopeia::KSParticle &aFinalParticle, unsigned index, bool bStart );
            double calcOrbitPhase(double tX, double tY);
            double quadrantOrbitCorrection(double phase, double vx);




        private:
            std::deque<locust::Particle> fNewParticleHistory;
            double fPitchAngle;
            double fT0trapMin;
            int fNCrossings;
            FieldCalculator* fFieldCalculator;
            kl_interface_ptr_t fInterface;
            unsigned fSampleIndex;
            unsigned fStartingIndex;
    };


}

#endif

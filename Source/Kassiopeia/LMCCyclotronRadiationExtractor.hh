/*
 * LMCCyclotronRadiationExtractor.hh
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 * @brief It transfers particle track information from Kassiopiea to Locust for field calculation. 
 */

#ifndef LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_
#define LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_

#include "LMCParticle.hh"

#include "KSTrajectory.h"

#include "KField.h"
#include "KSStepModifier.h"
#include "KSComponentTemplate.h"

#include "KToolbox.h"


namespace locust
{

    class CyclotronRadiationExtractor :
            public Kassiopeia::KSComponentTemplate< CyclotronRadiationExtractor, Kassiopeia::KSStepModifier >
    {
        public:
            CyclotronRadiationExtractor();
            CyclotronRadiationExtractor( const CyclotronRadiationExtractor& aOrig );
            virtual ~CyclotronRadiationExtractor();

            CyclotronRadiationExtractor* Clone() const;

        public:
            double e_prestep=0.;
            double e_poststep=0.;
            double t_prestep=0.;
            double de_step=0.;
            double dt_step=0.;
            double de = 0.;
            double dt = 0.;

            bool ExecutePreStepModification( Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticleQueue& aQueue );
            bool ExecutePostStepModification( Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle, Kassiopeia::KSParticleQueue& aQueue );

            void SetTrajectory( Kassiopeia::KSTrajectory* aTrajectory );
            void SetP8Phase( int P8Phase );


        private:
            void InitializeComponent();
            void DeinitializeComponent();
            double GetGroupVelocityTE11(Kassiopeia::KSParticle& aFinalParticle);
            double GetGroupVelocityTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetGroupVelocityTE01(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase2(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase1(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTE11(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTE01(Kassiopeia::KSParticle& aFinalParticle);
            double GetTM01FieldAfterBounces(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);
            double GetTE11FieldAfterOneBounce(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);
            double GetTE01FieldAfterOneBounce(Kassiopeia::KSParticle& anInitialParticle, Kassiopeia::KSParticle& aFinalParticle);
            locust::Particle ExtractKassiopeiaParticle( Kassiopeia::KSParticle &aFinalParticle);
            Kassiopeia::KSTrajectory* fTrajectory;
            int fP8Phase; // 1, 2, 3, or 4.
            std::deque<locust::Particle> fNewParticleHistory;



        protected:
            virtual void PullDeupdateComponent();
            virtual void PushDeupdateComponent();

    };

} /* namespace locust */

#endif /* LOCUST_LMCCYCLOTRONRADIATIONEXTRACTOR_HH_ */

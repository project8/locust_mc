#ifndef LMCFIELDCALCULATOR_HH_
#define LMCFIELDCALCULATOR_HH_


#include "KSSpaceInteraction.h"
#include "KSList.h"

#include "KSStep.h"
#include "KSTrajectory.h"

#include "KMathBracketingSolver.h"

#include "LMCConst.hh"
#include "LMCThreeVector.hh"
#include <vector>
using std::vector;





using katrin::KMathBracketingSolver;

namespace locust
{

    class FieldCalculator :
        public Kassiopeia::KSComponentTemplate< FieldCalculator, Kassiopeia::KSSpaceInteraction >
    {
        public:
            FieldCalculator();
            FieldCalculator( const FieldCalculator& aCopy );
            FieldCalculator* Clone() const;
            ~FieldCalculator();

            //*****************
            //space interaction
            //*****************



            double GetGroupVelocityTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetGroupVelocityTE01(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase2(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase1(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTE01(Kassiopeia::KSParticle& aFinalParticle);
            double GetTM01FieldWithTerminator(Kassiopeia::KSParticle& aFinalParticle);
            double GetTE01FieldAfterOneBounce(Kassiopeia::KSParticle& aFinalParticle);
      double GetRetardedTE01FieldAfterOneBounce(Kassiopeia::KSParticle& aFinalParticle);

            int FindNode(double tNew) const;
      double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition, double GroupVelocity);





        public:
      void CalculateInteraction( const Kassiopeia::KSTrajectory& aTrajectory, const Kassiopeia::KSParticle& aTrajectoryInitialParticle, const Kassiopeia::KSParticle& aTrajectoryFinalParticle, const KThreeVector& aTrajectoryCenter, const double& aTrajectoryRadius, const double& aTrajectoryTimeStep, Kassiopeia::KSParticle& anInteractionParticle, double& aTimeStep, bool& aFlag );
            void ExecuteInteraction( const Kassiopeia::KSParticle& anInteractionParticle, Kassiopeia::KSParticle& aFinalParticle, Kassiopeia::KSParticleQueue& aSecondaries ) const;

            //***********
            //composition
            //***********

        public:
            void AddSpaceInteraction( Kassiopeia::KSSpaceInteraction* anInteraction );
            void RemoveSpaceInteraction( Kassiopeia::KSSpaceInteraction* anInteraction );

        private:
            Kassiopeia::KSList< Kassiopeia::KSSpaceInteraction > fSpaceInteractions;
            Kassiopeia::KSSpaceInteraction* fSpaceInteraction;

            //******
            //action
            //******

        public:
            void SetStep( Kassiopeia::KSStep* anStep );
            void SetTrajectory( Kassiopeia::KSTrajectory* aTrajectory );

            void CalculateInteraction();
            void ExecuteInteraction();

            virtual void PushUpdateComponent();
            virtual void PushDeupdateComponent();

        private:
            Kassiopeia::KSStep* fStep;
            const Kassiopeia::KSParticle* fTerminatorParticle;
            const Kassiopeia::KSParticle* fTrajectoryParticle;
            Kassiopeia::KSParticle* fInteractionParticle;
            Kassiopeia::KSParticle* fFinalParticle;
            Kassiopeia::KSParticleQueue* fParticleQueue;
            Kassiopeia::KSTrajectory* fTrajectory;
    };

}

#endif

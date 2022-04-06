#ifndef LMCFIELDCALCULATOR_HH_
#define LMCFIELDCALCULATOR_HH_

#include "LMCKassLocustInterface.hh"

#include "KSSpaceInteraction.h"
#include "KSList.h"

#include "KSStep.h"

namespace Kassiopeia
{
    class KSParticle;
}

namespace locust
{

    class FieldCalculator
    {
        public:
            FieldCalculator();
            FieldCalculator( const FieldCalculator& aCopy );
            FieldCalculator* Clone() const;
            ~FieldCalculator();

            double GetGroupVelocityTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetGroupVelocityTE10(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase2(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase1(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorCavity(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTE10(Kassiopeia::KSParticle& aFinalParticle);
            double GetTM01FieldWithTerminator(Kassiopeia::KSParticle& aFinalParticle);
            double GetTE10FieldAfterOneBounce(Kassiopeia::KSParticle& aFinalParticle);
            double GetCavityDotProductFactor(Kassiopeia::KSParticle& aFinalParticle, std::vector<double> aTE_E_normalized);
            std::vector<double> GetCavityNormalizedModeField(int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle);
            double GetCavityFIRSample(Kassiopeia::KSParticle& aFinalParticle, int nFilterBinsRequired, double dtFilter);
            double calcOrbitPhase(double vx, double vy);
            double calcTheta(double x, double y);
            double quadrantOrbitCorrection(double phase, double vx);
            double quadrantPositionCorrection(double phase, double x);

            kl_interface_ptr_t fInterface;
    };

}

#endif

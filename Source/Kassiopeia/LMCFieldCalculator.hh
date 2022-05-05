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
            std::vector<double> GetWaveguideNormalizedModeField(int l, int m, int n, std::vector<double> tKassParticleXP);
            std::vector<double> GetCavityNormalizedModeField(int l, int m, int n, std::vector<double> tLocation, bool TE, bool Electric);
            double GetCavityDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> anE_normalized, bool IntermediateFile);
            double GetWaveguideDotProductFactor(std::vector<double> tKassParticleXP, std::vector<double> aTE_E_normalized, bool IntermediateFile);
            double GetCavityFIRSample(std::vector<double> tKassParticleXP, bool BypassTF);

            kl_interface_ptr_t fInterface;
    };

}

#endif

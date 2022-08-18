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
            double GetCouplingFactorTE011Cavity(Kassiopeia::KSParticle& aFinalParticle);
            double GetTM01FieldWithTerminator(Kassiopeia::KSParticle& aFinalParticle);
            double GetTE10FieldAfterOneBounce(Kassiopeia::KSParticle& aFinalParticle);
            double GetTE011FieldCavity(Kassiopeia::KSParticle& aFinalParticle);
            std::pair<double,double> GetCavityFIRSample(std::vector<double> tKassParticleXP, bool BypassTF);
            void SetNFilterBinsRequired( int aNumberOfBins );
            int GetNFilterBinsRequired();
            void SetFilterSize( int aFilterSize );

            kl_interface_ptr_t fInterface;

        private:
            std::deque<double> fFIRBuffer;
            std::deque<double> fFrequencyBuffer;
            int fNFilterBinsRequired;
    };

}

#endif

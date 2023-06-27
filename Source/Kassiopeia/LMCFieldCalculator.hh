#ifndef LMCFIELDCALCULATOR_HH_
#define LMCFIELDCALCULATOR_HH_

#include "LMCKassLocustInterface.hh"
#include "LMCEquivalentCircuit.hh" // : LMCAnalyticResponseFunction
#include "LMCDampedHarmonicOscillator.hh" // : LMCAnalyticResponseFunction
#include "LMCFIRFileHandler.hh"
#include "LMCTFFileHandler.hh"


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
            bool Configure( const scarab::param_node& aParam );
            bool ConfigureByInterface();


            bool ModeSelect(int l, int m, int n, bool eGun, bool bNormCheck, bool bTE);
            double GetGroupVelocityTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetGroupVelocityTE10(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase2(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase1(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorCavity(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTE10(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTXlmnCavity(int l, int m, int n, bool bTE, Kassiopeia::KSParticle& aFinalParticle);
            double GetTM01FieldWithTerminator(Kassiopeia::KSParticle& aFinalParticle);
            double GetTE10FieldAfterOneBounce(Kassiopeia::KSParticle& aFinalParticle);
            double GetTXlmnFieldCavity(int l, int m, int n, bool bTE, Kassiopeia::KSParticle& aFinalParticle);
            std::pair<double,double> GetCavityFIRSample(std::vector<double> tKassParticleXP, bool BypassTF);
            void SetNFilterBinsRequired( double dt );
            int GetNFilterBinsRequired();
            void SetFilterSize( int aFilterSize );
            int GetFilterSize();


            kl_interface_ptr_t fInterface;

        private:
            TFReceiverHandler* fTFReceiverHandler;
            AnalyticResponseFunction* fAnalyticResponseFunction;
            std::deque<double> fFIRBuffer;
            std::deque<double> fFrequencyBuffer;
            int fNFilterBinsRequired;
            bool fbMultiMode;
    };

}

#endif

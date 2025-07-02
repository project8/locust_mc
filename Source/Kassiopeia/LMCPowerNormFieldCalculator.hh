#ifndef LMCPOWERNORMFIELDCALCULATOR_HH_
#define LMCPOWERNORMFIELDCALCULATOR_HH_

#include "LMCKassLocustInterface.hh"
#include "LMCFieldCalculator.hh"
#include "LMCEquivalentCircuit.hh" // : LMCAnalyticResponseFunction
#include "LMCDampedHarmonicOscillator.hh" // : LMCAnalyticResponseFunction
#include "LMCFIRFileHandler.hh"
#include "LMCTFFileHandler.hh"

#include "KThreeVector.hh"
using katrin::KThreeVector;


#include "KSSpaceInteraction.h"
#include "KSList.h"

#include "KSStep.h"

namespace Kassiopeia
{
    class KSParticle;
}

namespace locust
{

    class PowerNormFieldCalculator : public FieldCalculator
    {
        public:
            PowerNormFieldCalculator();
            PowerNormFieldCalculator( const PowerNormFieldCalculator& aCopy );
            PowerNormFieldCalculator* Clone() const;
            ~PowerNormFieldCalculator();
            bool Configure( const scarab::param_node& aParam );
            bool ConfigureByInterface();


//            double GetGroupVelocityTM01(Kassiopeia::KSParticle& aFinalParticle);
//            double GetGroupVelocityTE10(Kassiopeia::KSParticle& aFinalParticle);
//            double GetDampingFactorPhase2(Kassiopeia::KSParticle& aFinalParticle);
//            double GetDampingFactorPhase1(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorCavity(Kassiopeia::KSParticle& aFinalParticle);
//            double GetCouplingFactorTM01(Kassiopeia::KSParticle& aFinalParticle);
//            double GetCouplingFactorTE10(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTXlmnCavity(int l, int m, int n, bool bTE, Kassiopeia::KSParticle& aFinalParticle);
//            double GetTM01FieldWithTerminator(Kassiopeia::KSParticle& aFinalParticle);
//            double GetTE10FieldAfterOneBounce(Kassiopeia::KSParticle& aFinalParticle);
            double GetTXlmnFieldCavity(int l, int m, int n, bool bTE, Kassiopeia::KSParticle& aFinalParticle);

            void SetCavityFIRSample(int bTE, int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle, bool BypassTF);
            std::pair<double,double> GetCavityFIRSample(int bTE, int l, int m, int n);

            KThreeVector InterpolateVelocity(double dt, double tCyclotronFrequency, double tCyclotronRadius, KThreeVector tVelocityParallel, KThreeVector fMagneticField, KThreeVector tAlpha, KThreeVector tBeta);
            KThreeVector InterpolatePosition(double dt, double tCyclotronFrequency, double tCyclotronRadius, KThreeVector tVelocityParallel, KThreeVector tMagneticField, KThreeVector tGuidingCenterPosition, KThreeVector tAlpha, KThreeVector tBeta);
            double calcTheta(double x, double y);

            void SetNFilterBinsRequired( double dt );
//            int GetNFilterBinsRequired();
            // void SetFilterSize( int aFilterSize );
            int GetFilterSize();


            kl_interface_ptr_t fInterface;

        private:
            double calcOrbitPhase(double vx, double vy);
            double quadrantOrbitCorrection(double phase, double vx);
            TFReceiverHandler* fTFReceiverHandler;
            AnalyticResponseFunction* fAnalyticResponseFunction;
            std::deque<double> fJdotEBuffer;
            int fNFilterBinsRequired;
            std::vector<std::vector<int>> fModeSet;
            double fTime;

    };

}

#endif

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

            bool ReconvertAnalyticGFtoFIR(int l, int m, int n);
            double GetGroupVelocityTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetGroupVelocityTE10(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase2(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorPhase1(Kassiopeia::KSParticle& aFinalParticle);
            double GetDampingFactorCavity(int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTM01(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorTE10(Kassiopeia::KSParticle& aFinalParticle);
            double GetCouplingFactorCavity(int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle);
            double GetTM01FieldWithTerminator(Kassiopeia::KSParticle& aFinalParticle);
            double GetTE10FieldAfterOneBounce(Kassiopeia::KSParticle& aFinalParticle);
            double GetFieldCavity(int l, int m, int n, Kassiopeia::KSParticle& aFinalParticle);
            std::pair<double,double> GetCavityFIRSample(int l, int m, int n, std::vector<double> tKassParticleXP, bool BypassTF);
            void SetNFilterBinsRequired( double dt );
            int GetNFilterBinsRequired();
            void SetFilterSize( int aFilterSize );
            int GetFilterSize();
            void SetNFilterBinsRequiredArray(int l, int m, int n, double dt );
            int GetNFilterBinsRequiredArray(int l, int m, int n);
            void SetFilterSizeArray(int l, int m, int n, int aFilterSize );
            int GetFilterSizeArray(int l, int m, int n);

            kl_interface_ptr_t fInterface;

        private:
            TFReceiverHandler* fTFReceiverHandler;
            AnalyticResponseFunction* fAnalyticResponseFunction;
            std::deque<double> fFIRBuffer;
            std::deque<double> fFrequencyBuffer;
            std::vector< std::vector< std::vector< std::deque<double> >>> fFIRBufferArray;
            std::vector< std::vector< std::vector< std::deque<double> >>> fFrequencyBufferArray;
            int fNFilterBinsRequired;
	    std::vector< std::vector< std::vector< int >>> fNFilterBinsRequiredArray;
    };

}

#endif

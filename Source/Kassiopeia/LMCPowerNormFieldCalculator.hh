#ifndef LMCPOWERNORMFIELDCALCULATOR_HH_
#define LMCPOWERNORMFIELDCALCULATOR_HH_

#include "LMCKassLocustInterface.hh"
#include "LMCFieldCalculator.hh"
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

    class PowerNormFieldCalculator : public FieldCalculator
    {
        public:
            PowerNormFieldCalculator();
            PowerNormFieldCalculator( const PowerNormFieldCalculator& aCopy );
            PowerNormFieldCalculator* Clone() const;
            ~PowerNormFieldCalculator();

            kl_interface_ptr_t fInterface;


        private:
            double calcOrbitPhase(double vx, double vy);
            double quadrantOrbitCorrection(double phase, double vx);
            TFReceiverHandler* fTFReceiverHandler;
            AnalyticResponseFunction* fAnalyticResponseFunction;
            std::deque<double> fFIRBuffer;
            std::deque<double> fFrequencyBuffer;
            int fNFilterBinsRequired;
            std::vector<std::vector<int>> fModeSet;

    };

}

#endif

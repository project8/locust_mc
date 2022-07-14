/*
 * LMCDampedHarmonicOscillator.hh
 *
 *  Created on: Jul 6, 2022
 *      Author: pslocum
 *
This class calculates an analytic Green's function from a model of a damped harmonic
oscillator.  To use the model in the context of the LMCCavitySignalGenerator, the
following default syntax can be used:

LocustSim -c config/LocustCavityTemplate.json

To print out the resulting FIR to a text file "output/FIR.txt", use this parameter:

"cavity-signal.print-fir-debug"=true
 *
 *
 *
 */

#ifndef LMCDAMPEDHARMONICOSCILLATOR_HH_
#define LMCDAMPEDHARMONICOSCILLATOR_HH_
#include "LMCAnalyticResponseFunction.hh"
#include "LMCConst.hh"
#include "param.hh"
#include "LMCException.hh"
#include <math.h>

namespace locust
{
    class DampedHarmonicOscillator : public AnalyticResponseFunction
    {

        public:
    		DampedHarmonicOscillator();
    		virtual ~DampedHarmonicOscillator();
    		virtual bool Configure( const scarab::param_node& aNode );
            virtual bool GenerateGreensFunction();
            bool Initialize();
            double GreensFunction(double t);
            double ExpDecayTerm(double t);


        private:
            double fCavityFrequency; // Hz
            double fCavityOmega; // radians/s
            double fCavityQ;
            double fCavityDampingFactor;
            double fBFactor; // harmonic oscillator parameter
            double fCavityOmegaPrime;  // damped resonant frequency
            int fMaxNBins;
            double fTimeResolution;
            double fThresholdFactor;

    };


} /* namespace locust */

#endif

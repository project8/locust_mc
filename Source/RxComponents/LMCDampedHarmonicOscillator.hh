/*
 * LMCDampedHarmonicOscillator.hh
 *
 *  Created on: Jul 6, 2022
 *      Author: pslocum
 */

#ifndef LMCEQUIVALENTCIRCUIT_HH_
#define LMCEQUIVALENTCIRCUIT_HH_
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
            virtual void GenerateFIR();
            bool Initialize();
            double GreensFunction(double t);
            double ExpDecayTerm(double t);
            std::vector<std::pair<double,double>> gfArray;


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

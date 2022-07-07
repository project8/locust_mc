/*
 * LMCDampedHarmonicOscillator.hh
 *
 *  Created on: Jul 6, 2022
 *      Author: pslocum
 */

#ifndef LMCEQUIVALENTCIRCUIT_HH_
#define LMCEQUIVALENTCIRCUIT_HH_
#include "LMCAnalyticResponseFunction.hh"
#include "param.hh"
#include "LMCException.hh"

namespace locust
{
    class DampedHarmonicOscillator : public AnalyticResponseFunction
    {

        public:
    		DampedHarmonicOscillator();
    		virtual ~DampedHarmonicOscillator();
    		virtual bool Configure( const scarab::param_node& aNode );
    };


} /* namespace locust */

#endif

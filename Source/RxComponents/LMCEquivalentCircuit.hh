/*
 * LMCEquivalentCircuit.hh
 *
 *  Created on: Nov 12, 2021
 *      Author: jgaison

This class calculates an analytic transfer function from a model of an equivalent
resonant RLC circuit.  The default resonant frequency is 1 GHz with Q=1000.  To use
the model in the context of the LMCCavitySignalGenerator, the following parameters
are needed in the configuration, either in the json file or on the command line:

"cavity-signal.equivalent-circuit"=true,
"cavity-signal.shift-n-bins": 6000,
"cavity-signal.zero-padding-size": 1000

To print out the resulting FIR to a text file "output/FIR.txt", use this parameter:

"cavity-signal.print-fir-debug"=true

 */

#ifndef LMCEQUIVALENTCIRCUIT_HH_
#define LMCEQUIVALENTCIRCUIT_HH_
#include "LMCAnalyticResponseFunction.hh"
#include "param.hh"
#include "LMCException.hh"
#include "LMCConst.hh"
#include "LMCSignal.hh"
#include "LMCPatchAntenna.hh"
#include "LMCSlotAntenna.hh"
#include <vector>


namespace locust
{
    class EquivalentCircuit : public AnalyticResponseFunction
    {

        public:
    		EquivalentCircuit();
    		virtual ~EquivalentCircuit();
    		virtual bool Configure( const scarab::param_node& aNode );
    		virtual bool GenerateTransferFunction();

        private:
    		int nbins;
    		double fEquivalentR;
    		double fEquivalentL;
    		double fEquivalentC;
    		int fTFBins;
    		double fFreqRangeCenter;
};


} /* namespace locust */

#endif

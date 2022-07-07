/*
 * LMCEquivalentCircuit.hh
 *
 *  Created on: Nov 12, 2021
 *      Author: jgaison
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
#include "LMCComplexFFT.hh"
#include <vector>


namespace locust
{
    class EquivalentCircuit : public AnalyticResponseFunction
    {

        public:
    		EquivalentCircuit();
    		virtual ~EquivalentCircuit();
    		virtual bool Configure( const scarab::param_node& aNode );
    		void GenerateTransferFunction();
    		std::vector<std::complex<double>> tfArray;
    		double initialFreq;
    		bool fGeneratingTF;

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

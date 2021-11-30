/*
 * LMCEquivalentCircuit.hh
 *
 *  Created on: Nov 12, 2021
 *      Author: jgaison
 */

#ifndef LMCEQUIVALENTCIRCUIT_HH_
#define LMCEQUIVALENTCIRCUIT_HH_
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
 /*!
 @class LMCEquivalentCicuit
 @author J. Gaison
 @brief Base class to implement cavity parameterization
 @details
 Available configuration options:
 R, L, and C of effective circuit for cavity
 */
    class EquivalentCircuit
    {

        public:
            EquivalentCircuit();
            virtual ~EquivalentCircuit();
	    void GenerateTransferFunction(double R, double L, double C);
	    std::vector<std::complex<double>> tfArray;
	    double initialFreq;

        private:
	    int nbins;
};


} /* namespace locust */

#endif

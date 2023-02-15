/*
 * LMCCavityUtility.hh
 *
 *  Created on: Feb 10, 2023
 *      Author: pslocum
 */

#ifndef LMCCAVITYUTILITY_HH_
#define LMCCAVITYUTILITY_HH_

#include "LMCSignal.hh"
#include "LMCTFFileHandler.hh"
#include "LMCFIRFileHandler.hh"
#include "LMCAnalyticResponseFunction.hh"
#include "LMCDampedHarmonicOscillator.hh"
#include "LMCUtility.hh"

namespace locust
{

    /*!
     @class CavityUtility
     @author P. L. Slocum Feb. 10 2023

     @brief Class to generate small cavity tests

     @details
     Operates in time space

     Configuration name: N/A

     Available configuration options: N/A

     */
    class CavityUtility : public Utility
    {
    public:

        CavityUtility();
        virtual ~CavityUtility();

        bool Configure();
        bool CheckCavityQ( double dhoTimeResolution, double dhoThresholdFactor, double dhoCavityFrequency, double dhoCavityQ);
    	void AddParam(std::string aString, double aValue);
    	std::deque<double> SignalToDeque(Signal* aSignal);
    	bool PopulateSignal(Signal* aSignal, int N0);
    	const scarab::param_node* GetParams();


    private:

        TFReceiverHandler* fTFReceiverHandler;
        AnalyticResponseFunction* fAnalyticResponseFunction;
        scarab::param_node* fparam_0;

        double fRF_frequency;
        double fFilterRate;


    };


} /* namespace locust */

#endif /* LMCCAVITYUTILITY_HH_ */

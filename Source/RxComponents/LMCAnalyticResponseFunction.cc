/*
 * LMCAnalyticResponseFunction.cc
 *
 *  Created on: Jul6, 2022
 *      Author: pslocum
 */

#include "LMCAnalyticResponseFunction.hh"
#include <iostream>
#include "logger.hh"

namespace locust
{

	LOGGER( lmclog, "AnalyticResponseFunction" );


    AnalyticResponseFunction::AnalyticResponseFunction():
	    fGeneratingTF( false ),
		fInitialFreq( 0. ),
		fTFarray(0 )
    {}
    AnalyticResponseFunction::~AnalyticResponseFunction() {}


    bool AnalyticResponseFunction::Configure( const scarab::param_node& aParam )
    {
        return true;
    }

    void AnalyticResponseFunction::SetGeneratingTF( bool aFlag )
    {
    	fGeneratingTF = aFlag;
    }
    bool AnalyticResponseFunction::GetGeneratingTF()
    {
    	return fGeneratingTF;
    }
    void AnalyticResponseFunction::SetInitialFreq( double aFreq )
    {
    	fInitialFreq = aFreq;
    }
    double AnalyticResponseFunction::GetInitialFreq()
    {
    	return fInitialFreq;
    }
    void AnalyticResponseFunction::SetTFarray( std::vector<std::complex<double>> aTFarray )
    {
    	fTFarray = aTFarray;
    }
    std::vector<std::complex<double>> AnalyticResponseFunction::GetTFarray()
    {
    	return fTFarray;
    }








} /* namespace locust */


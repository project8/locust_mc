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

    void AnalyticResponseFunction::SetGFarray( std::vector<std::pair<double,std::pair<double,double> > > aGFarray )
    {
    	if (fGFarray.size() > 0)
    	{
        	for (unsigned index=0; index<aGFarray.size(); index++)
        	{
        		fGFarray[index].second.first = aGFarray[index].second.first;
        		fGFarray[index].second.second = aGFarray[index].second.second;
        	}
    	}
    	else
    	{
    		for (unsigned index=0; index<aGFarray.size(); index++)
    		{
    			fGFarray.push_back(std::make_pair(aGFarray[index].first, aGFarray[index].second));
    		}
    	}
    }

    std::vector<std::pair<double,std::pair<double,double> > > AnalyticResponseFunction::GetGFarray()
    {
    	return fGFarray;
    }








} /* namespace locust */


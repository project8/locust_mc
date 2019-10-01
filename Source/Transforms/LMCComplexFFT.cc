/*
 * LMCComplexFFT.cc
 *
 *  Created on: Sep 30, 2019
 *      Author: P. T. Surukuchi
 */

#include "LMCComplexFFT.hh"

namespace locust
{
    ComplexFFT::ComplexFFT()
    {
    }

    ComplexFFT::~ComplexFFT()
    {
    }

    bool ComplexFFT::Configure(const scarab::param_node& aParam)
    {
    	if( aParam.has("Transform-flag"))
        {
    		fbufferMargin=aParam["Transform-flag"]().as_bool();
        }
       	if(aParam.has("use-wisdom"))
        {
       		fbufferSize=aParam["use-wisdom"]().as_bool();
        }
        if(aParam.has("wisdom-filename"))
        {
            fbufferSize=aParam["wisdom-filename"]().as_string();
        }
    	return true;
    }
} /* namespace locust */

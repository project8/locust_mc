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


    AnalyticResponseFunction::AnalyticResponseFunction()
    {}
    AnalyticResponseFunction::~AnalyticResponseFunction() {}


    bool AnalyticResponseFunction::Configure( const scarab::param_node& aParam )
    {
        return true;
    }



} /* namespace locust */


/*
 * LMCField.cc
 *
 *  Created on: Jun 4, 2021
 *      Author: pslocum
 */

#include "LMCField.hh"


namespace locust
{
    LOGGER( lmclog, "Field" );
    Field::Field() {}
    Field::~Field() {}


    std::vector<double> Field::GetNormFactors()
    {
    	return fNormFactor;
    }

    void Field::SetNormFactors(std::vector<double> aNormFactor)
    {
    	fNormFactor = aNormFactor;
    }





} /* namespace locust */


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


    std::vector<double> Field::GetNormFactorsTE()
    {
    	return fNormFactorTE;
    }

    void Field::SetNormFactorsTE(std::vector<double> aNormFactor)
    {
    	fNormFactorTE = aNormFactor;
    }

    std::vector<double> Field::GetNormFactorsTM()
    {
    	return fNormFactorTM;
    }

    void Field::SetNormFactorsTM(std::vector<double> aNormFactor)
    {
    	fNormFactorTM = aNormFactor;
    }





} /* namespace locust */


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
    Field::Field():
    fCentralFrequency( 25.9e9 )
    {}
    Field::~Field() {}


    std::vector<std::vector<std::vector<double>>> Field::GetNormFactorsTE()
    {
    	return fModeNormFactorTE;
    }

    void Field::SetNormFactorsTE(std::vector<std::vector<std::vector<double>>> aNormFactor)
    {
    	fModeNormFactorTE = aNormFactor;
    }

    std::vector<std::vector<std::vector<double>>> Field::GetNormFactorsTM()
    {
    	return fModeNormFactorTM;
    }

    void Field::SetNormFactorsTM(std::vector<std::vector<std::vector<double>>> aNormFactor)
    {
    	fModeNormFactorTM = aNormFactor;
    }


    double Field::GetCentralFrequency()
    {
    	return fCentralFrequency;
    }

    void Field::SetCentralFrequency( double aCentralFrequency )
    {
    	fCentralFrequency = aCentralFrequency;
    }




} /* namespace locust */


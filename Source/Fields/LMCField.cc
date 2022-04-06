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
		fnPixels( 100 ),
    	fCentralFrequency(0.)
    {}
    Field::~Field() {}

    bool Field::Configure( const scarab::param_node& aParam )
    {

    	if( aParam.has( "central-frequency" ) )
        	{
            	fCentralFrequency= 2.*LMCConst::Pi()*aParam["central-frequency"]().as_double();
        	}

       	if( aParam.has( "n-pixels" ) )
            {
                SetNPixels(aParam["n-pixels"]().as_int());
            }

    	return true;

    }


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

    int Field::GetNPixels()
    {
    	return fnPixels;
    }

    void Field::SetNPixels( int aNumberOfPixels )
    {
    	fnPixels = aNumberOfPixels;
    }




} /* namespace locust */


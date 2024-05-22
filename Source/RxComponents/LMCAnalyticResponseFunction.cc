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
		fNModes( 2 ),
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
    void AnalyticResponseFunction::SetNModes( int aNumberOfModes )
    {
    	fNModes = aNumberOfModes;
    }
    int AnalyticResponseFunction::GetNModes()
    {
    	return fNModes;
    }
    void AnalyticResponseFunction::SetTFarray( std::vector<std::complex<double>> aTFarray )
    {
    	fTFarray = aTFarray;
    }
    std::vector<std::complex<double>> AnalyticResponseFunction::GetTFarray()
    {
    	return fTFarray;
    }

    void AnalyticResponseFunction::SetGFarray( std::vector <std::vector< std::vector< std::vector< std::vector<std::pair<double,std::pair<double,double> > > > > > > aGFarray )
    {
        fGFarray = aGFarray;
    }

    std::vector< std::vector<std::pair<double,std::pair<double,double> > > > AnalyticResponseFunction::GetGFarray( std::vector<std::vector<int>> aModeSet)
    {
        std::vector< std::vector<std::pair<double,std::pair<double,double>>>> anArrayOfGFArrays;
        for (int mu=0; mu<aModeSet.size(); mu++)
		{
		    bool bTE = aModeSet[mu][0];
		    int l = aModeSet[mu][1];
		    int m = aModeSet[mu][2];
		    int n = aModeSet[mu][3];

		    anArrayOfGFArrays.push_back( fGFarray[bTE][l][m][n] );
		}

        return anArrayOfGFArrays;
    }








} /* namespace locust */


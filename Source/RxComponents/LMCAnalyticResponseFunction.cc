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
    {  
	int fNModes = 2;
	fGFarray.resize(2); 
	for( int bTE = 0; bTE<2; bTE++)
	{
        	fGFarray[bTE].resize(fNModes);
        	for(int l=0; l<fNModes; l++)
        	{   
                	fGFarray[bTE][l].resize(fNModes);
                	for(int m=0; m<fNModes; m++)
                	{   
                        	fGFarray[bTE][l][m].resize(fNModes);
                	}
		}    
        }   
    } 
    AnalyticResponseFunction::~AnalyticResponseFunction() {}


    bool AnalyticResponseFunction::Configure( const scarab::param_node& aParam )
    {

	if( aParam.has( "n-modes" ) )
	{
		int fNModes = aParam["n-modes"]().as_int();
        	for( int bTE = 0; bTE<2; bTE++)
        	{   
                	fGFarray[bTE].resize(fNModes);
                	for(int l=0; l<fNModes; l++)
                	{   
                        	fGFarray[bTE][l].resize(fNModes);
                        	for(int m=0; m<fNModes; m++)
                        	{   
                                	fGFarray[bTE][l][m].resize(fNModes);
                        	}   
                	}    
        	}  
	}	
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
    void AnalyticResponseFunction::SetGFarray(int bTE, int l, int m, int n, std::vector<std::pair<double,std::pair<double,double> > > aGFarray )
    {
	if(fGFarray[bTE][l][m][n].size()!=aGFarray.size())
	{
		fGFarray[bTE][l][m][n].resize(aGFarray.size());
	}
	for (unsigned index=0; index<aGFarray.size(); index++)
        {   
        	fGFarray[bTE][l][m][n][index] = std::make_pair(aGFarray[index].first, aGFarray[index].second);
        }  
    }

    std::vector<std::pair<double,std::pair<double,double> > > AnalyticResponseFunction::GetGFarray(int bTE, int l, int m, int n)
    {
    	return fGFarray[bTE][l][m][n];
    }








} /* namespace locust */


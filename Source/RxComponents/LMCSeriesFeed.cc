/*
 * LMCSeriesFeed.cc
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#include "LMCSeriesFeed.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "SeriesFeed" );

    SeriesFeed::SeriesFeed()
    {
    }

    SeriesFeed::~SeriesFeed()
    {
    }

    bool SeriesFeed::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombinerParent::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from SeriesFeed subclass");
    	}

    	SetVoltageDampingFactors();
    	return true;
    }


	bool SeriesFeed::SetVoltageDampingFactors()
	{
		return true;

	}




} /* namespace locust */

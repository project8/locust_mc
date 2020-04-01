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

    	if( !PowerCombiner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from SeriesFeed subclass");
    		return false;
    	}

		SetJunctionLoss(0.87);
		SetPatchLoss(0.38);
		SetAmplifierLoss(0.66);
		SetEndPatchLoss(0.38);

    	SetVoltageDampingFactors();
    	return true;
    }


	bool SeriesFeed::SetVoltageDampingFactors()
	{
		for (unsigned z_index=0; z_index<GetNElementsPerStrip(); z_index++)
		{
			int njunctions = z_index;
			double aFactor = GetPatchLoss()*GetAmplifierLoss()*pow(GetJunctionLoss(), njunctions);
			SetDampingFactor(z_index, aFactor);
		}
		return true;
	}




} /* namespace locust */

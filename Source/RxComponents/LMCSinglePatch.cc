/*
 * LMCSinglePatch.cc
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#include "LMCSinglePatch.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "SinglePatch" );

    SinglePatch::SinglePatch()
    {
    }

    SinglePatch::~SinglePatch()
    {
    }

    bool SinglePatch::Configure( const scarab::param_node& aParam )
    {

    	if( !PowerCombiner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring PowerCombiner class from SinglePatch subclass");
    	}
    	if (GetNElementsPerStrip() != 1 )
    	{
    		LERROR(lmclog,"single-patch power combining config expects one patch per strip.");
    		return false;
    	}

    	SetAmplifierLoss(0.7071);  // there is a single impedance transformation to the amplifier.
    	SetVoltageDampingFactors();
    	return true;
    }


    bool SinglePatch::IsSinglePatch()
    {
    	return true;
    }


	bool SinglePatch::SetVoltageDampingFactors()
	{

		for (unsigned z_index=0; z_index<GetNElementsPerStrip(); z_index++)
		{
			double aFactor = GetAmplifierLoss();
			SetDampingFactor(z_index, aFactor);
		}

		return true;

	}




} /* namespace locust */

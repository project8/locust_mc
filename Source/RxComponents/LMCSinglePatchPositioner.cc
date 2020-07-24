/*
 * LMCSinglePatchPositioner.cc
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#include "LMCSinglePatchPositioner.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "SinglePatchPositioner" );

    SinglePatchPositioner::SinglePatchPositioner()
    {
    }

    SinglePatchPositioner::~SinglePatchPositioner()
    {
    }

    bool SinglePatchPositioner::Configure( const scarab::param_node& aParam )
    {

    	if( !AntennaElementPositioner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring AntennaElementPositioner class from SinglePatchPositioner subclass");
    	}

    	return true;
    }

    double SinglePatchPositioner::GetPositionZ(double zShiftArray, int channelIndex, int nChannels,
    		int nSubarrays, int nReceivers, double elementSpacingZ, int receiverIndex)
    {
    	double zPosition = 0.;
    	return zPosition;
    }




} /* namespace locust */

/*
 * LMCPlanarArrayPositioner.cc
 *
 *  Created on: July 23, 2020
 *      Author: atelles and pslocum
 */

#include "LMCPlanarArrayPositioner.hh"
using std::string;


namespace locust
{

	LOGGER( lmclog, "PlanarArrayPositioner" );

    PlanarArrayPositioner::PlanarArrayPositioner()
    {
    }

    PlanarArrayPositioner::~PlanarArrayPositioner()
    {
    }

    bool PlanarArrayPositioner::Configure( const scarab::param_node& aParam )
    {

    	if( !AntennaElementPositioner::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring AntennaElementPositioner class from PlanarArrayPositioner subclass");
    	}

    	return true;
    }




} /* namespace locust */

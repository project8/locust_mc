/*
 * LMCAntennaElementPositioner.cc
 *
 *  Created on: Jul 22, 2020
 *      Author: pslocum
 */

#include "LMCAntennaElementPositioner.hh"
#include <iostream>
#include "logger.hh"

namespace locust
{

	LOGGER( lmclog, "AntennaElementPositioner" );


    AntennaElementPositioner::AntennaElementPositioner()
    {}
    AntennaElementPositioner::~AntennaElementPositioner() {}


    bool AntennaElementPositioner::Configure( const scarab::param_node& aParam )
    {
        return true;
    }




} /* namespace locust */


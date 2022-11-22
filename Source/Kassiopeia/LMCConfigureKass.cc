/*
 * LMCConfigureKass.cc
 *
 *  Created on: Oct 20, 2022
 *      Author: pslocum
 */

#include "LMCConfigureKass.hh"
#include "logger.hh"


namespace locust
{

    ConfigureKass::ConfigureKass() :
    		fParam(0)
    {
    }

    ConfigureKass::~ConfigureKass()
    {
    }

    void ConfigureKass::SetParameters( const scarab::param_node& aParam )
    {
    	fParam = &aParam;
    }


    const scarab::param_node* ConfigureKass::GetParameters()
    {
    	return fParam;
    }


} /* namespace locust */

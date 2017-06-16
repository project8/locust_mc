/*
 * LMCDefaultConfig.cc
 *
 *  Created on: Nov 4, 2013
 *      Author: nsoblath
 */

#include "LMCDefaultConfig.hh"

#include<string>
using std::string;

namespace locust
{

    DefaultConfig::DefaultConfig()
    {
        // default client configuration

        scarab::param_value tValue;

        add( "", tValue <<  );

        add( "", tValue <<  );

        //add( "port", tValue << 98342 );

        //add( "host", tValue << "localhost" );

        //add( "client-port", tValue << 98343 );

        //add( "client-host", tValue << "localhost" );

        //add( "file", tValue << "mantis_client_out.egg" );

        //add( "rate", tValue << 250.0 );

        //add( "duration", tValue << 1000 );

        //add( "mode", tValue << 0 );

        //add( "file-writer", tValue << "server" );
    }

    DefaultConfig::~DefaultConfig()
    {
    }

} /* namespace locust */

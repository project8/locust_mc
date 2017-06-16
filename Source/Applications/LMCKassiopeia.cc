/*
 * LMCKassiopeia.cc
 *
 *  Created on: Mar 10, 2016
 *      Author: nsoblath
 */

#include "LMCRunKassiopeia.hh"

#include "KCommandLineTokenizer.hh"

#include "KSMainMessage.h"


using namespace Kassiopeia;
using namespace katrin;
using namespace locust;

int main( int argc, char** argv )
{
    if( argc == 1 )
    {
        mainmsg( eNormal ) << "usage: ./Kassiopeia <config_file_one.xml> [<config_file_one.xml> <...>] [ -r variable1=value1 variable2=value ... ]" << eom;
        exit( -1 );
    }

    KCommandLineTokenizer tCommandLine;
    tCommandLine.ProcessCommandLine( argc, argv );

    RunKassiopeia tRunKass;
    tRunKass.SetVariableMap( tCommandLine.GetVariables() );

    int tReturn = tRunKass.Run( tCommandLine.GetFiles() );

    mainmsg( eNormal ) << "...finished" << eom;

    return tReturn;
}

/*
 * LocustMC.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#include "LMCConfigurator.hh"
//#include "LMCDefaultConfig.hh"
#include "LMCException.hh"
#include "LMCGeneratorToolbox.hh"
#include "LMCLogger.hh"
#include "LMCSimulationController.hh"

using namespace locust;

LMCLOGGER( lmclog, "LocustSim" );

int main( int argc, char** argv )
{
    LMCDEBUG( lmclog, "Welcome to Locust_MC\n\n" <<
            " (                                     *            \n" <<
            " )\\ )                         )      (  `      (    \n" <<
            "(()/(              (       ( /(      )\\))(     )\\   \n" <<
            " /(_))  (    (    ))\\  (   )\\())    ((_)()\\  (((_)  \n" <<
            "(_))    )\\   )\\  /((_) )\\ (_))/     (_()((_) )\\___  \n" <<
            "| |    ((_) ((_)(_))( ((_)| |_      |  \\/  |((/ __| \n" <<
            "| |__ / _ \\/ _| | || |(_-<|  _|     | |\\/| | | (__  \n" <<
            "|____|\\___/\\__|  \\_,_|/__/ \\__|_____|_|  |_|  \\___| \n" <<
            "                              |_____|               \n");

    Configurator* configurator = NULL;
    try
    {
        configurator = new Configurator( argc, argv );
    }
    catch( Exception& e )
    {
        LMCERROR( lmclog, "unable to configure LocustMC: " << e.what() );
        return -1;
    }

    LMCINFO( lmclog, "Setting up generator toolbox" );
    GeneratorToolbox toolbox;
    toolbox.Configure( configurator->config() );

    LMCINFO( lmclog, "Setting up simulation controller" );
    SimulationController controller;
    controller.Configure( configurator->config()->node_at( "simulation" ) );
    controller.SetFirstGenerator( toolbox.GetFirstGenerator() );

    LMCINFO( lmclog, "Beginning simulation run" );
    controller.Run();

    delete configurator;

    return 0;
}

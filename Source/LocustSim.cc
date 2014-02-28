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
            "\t\t (                                     *            \n" <<
            "\t\t )\\ )                         )      (  `      (    \n" <<
            "\t\t(()/(              (       ( /(      )\\))(     )\\   \n" <<
            "\t\t /(_))  (    (    ))\\  (   )\\())    ((_)()\\  (((_)  \n" <<
            "\t\t(_))    )\\   )\\  /((_) )\\ (_))/     (_()((_) )\\___  \n" <<
            "\t\t| |    ((_) ((_)(_))( ((_)| |_      |  \\/  |((/ __| \n" <<
            "\t\t| |__ / _ \\/ _| | || |(_-<|  _|     | |\\/| | | (__  \n" <<
            "\t\t|____|\\___/\\__|  \\_,_|/__/ \\__|_____|_|  |_|  \\___| \n" <<
            "\t\t                              |_____|               \n");

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
    toolbox.Configure( configurator->Config() );

    LMCINFO( lmclog, "Setting up simulation controller" );
    SimulationController controller;
    controller.SetFirstGenerator( toolbox.GetFirstGenerator() );
    controller.Configure( configurator->Config()->NodeAt( "simulation" ) );

    LMCINFO( lmclog, "Beginning simulation run" );
    controller.Run();

    delete configurator;

    return 0;
}

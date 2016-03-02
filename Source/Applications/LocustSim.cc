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
#include "LMCSimulationController.hh"

#include "logger.hh"

using namespace locust;

LOGGER( lmclog, "LocustSim" );

int main( int argc, char** argv )
{
    try
    {
        DEBUG( lmclog, "Welcome to Locust_MC\n\n" <<
                "\t\t (                                     *            \n" <<
                "\t\t )\\ )                         )      (  `      (    \n" <<
                "\t\t(()/(              (       ( /(      )\\))(     )\\   \n" <<
                "\t\t /(_))  (    (    ))\\  (   )\\())    ((_)()\\  (((_)  \n" <<
                "\t\t(_))    )\\   )\\  /((_) )\\ (_))/     (_()((_) )\\___  \n" <<
                "\t\t| |    ((_) ((_)(_))( ((_)| |_      |  \\/  |((/ __| \n" <<
                "\t\t| |__ / _ \\/ _| | || |(_-<|  _|     | |\\/| | | (__  \n" <<
                "\t\t|____|\\___/\\__|  \\_,_|/__/ \\__|_____|_|  |_|  \\___| \n" <<
                "\t\t                              |_____|               \n");

        Configurator configurator( argc, argv );

        INFO( lmclog, "Setting up generator toolbox" );
        GeneratorToolbox toolbox;
        if( ! toolbox.Configure( configurator.Config() ) )
        {
            ERROR( lmclog, "Unable to configure the generator toolbox" );
            return -1;
        }

        INFO( lmclog, "Setting up simulation controller" );
        SimulationController controller;
        controller.SetFirstGenerator( toolbox.GetFirstGenerator() );
        if( ! controller.Configure( configurator.Config()->NodeAt( "simulation" ) ) )
        {
            ERROR( lmclog, "Unable to configure the simulation controller" );
            return -1;
        }

        INFO( lmclog, "Preparing for run" );
        if( ! controller.Prepare() )
        {
            ERROR( lmclog, "Unable to prepare for the run" );
            return -1;
        }

        INFO( lmclog, "Beginning simulation run" );
        controller.Run();

        INFO( lmclog, "Run complete; finalizing" );
        controller.Finalize();
    }
    catch( Exception& e )
    {
        ERROR( lmclog, "Locust exception caught: " << e.what() );
        return -1;
    }
    catch( std::exception& e )
    {
        ERROR( lmclog, "Exception caught: " << e.what() );
        return -1;
    }

    return 0;
}

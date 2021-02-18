/*
 * LocustMC.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#include "configurator.hh"
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
        LPROG( lmclog, "Welcome to Locust_MC\n\n" <<
                "\t\t (                                     *            \n" <<
                "\t\t )\\ )                         )      (  `      (    \n" <<
                "\t\t(()/(              (       ( /(      )\\))(     )\\   \n" <<
                "\t\t /(_))  (    (    ))\\  (   )\\())    ((_)()\\  (((_)  \n" <<
                "\t\t(_))    )\\   )\\  /((_) )\\ (_))/     (_()((_) )\\___  \n" <<
                "\t\t| |    ((_) ((_)(_))( ((_)| |_      |  \\/  |((/ __| \n" <<
                "\t\t| |__ / _ \\/ _| | || |(_-<|  _|     | |\\/| | | (__  \n" <<
                "\t\t|____|\\___/\\__|  \\_,_|/__/ \\__|_____|_|  |_|  \\___| \n" <<
                "\t\t                              |_____|               \n");

        scarab::configurator configurator( argc, argv );

        LPROG( lmclog, "Setting up generator toolbox" );
        GeneratorToolbox toolbox;
        if( ! toolbox.Configure( configurator.config() ) )
        {
            LERROR( lmclog, "Unable to configure the generator toolbox" );
            return -1;
        }

        LPROG( lmclog, "Setting up simulation controller" );
        SimulationController controller;
        controller.SetFirstGenerator( toolbox.GetFirstGenerator() );
        if( ! controller.Configure( configurator.config()[ "simulation"].as_node()  ) )
        {
            LERROR( lmclog, "Unable to configure the simulation controller" );
            return -1;
        }

        LPROG( lmclog, "Preparing for run" );
        if( ! controller.Prepare() )
        {
            LERROR( lmclog, "Unable to prepare for the run" );
            return -1;
        }

        LPROG( lmclog, "Beginning simulation run" );
        controller.Run();

        LPROG( lmclog, "Run complete; finalizing" );
        controller.Finalize();
    }
    catch( Exception& e )
    {
        LERROR( lmclog, "Locust exception caught: " << e.what() );
        return -1;
    }
    catch( std::exception& e )
    {
        LERROR( lmclog, "Exception caught: " << e.what() );
        return -1;
    }
    catch( int &e )
    {
        LERROR( lmclog, "Exception caught: " << e );
        return e;
    }

    return 0;
}

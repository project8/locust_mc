/*
 * LocustMC.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#include "../Core/LMCConfigurator.hh"
//#include "LMCDefaultConfig.hh"
#include "../Core/LMCException.hh"
#include "../Core/LMCGeneratorToolbox.hh"
#include "../Core/LMCLogger.hh"
#include "../Core/LMCSimulationController.hh"

using namespace locust;

LMCLOGGER( lmclog, "LocustSim" );

int main( int argc, char** argv )
{
    try
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

        Configurator configurator( argc, argv );

        LMCINFO( lmclog, "Setting up generator toolbox" );
        GeneratorToolbox toolbox;
        if( ! toolbox.Configure( configurator.Config() ) )
        {
            LMCERROR( lmclog, "Unable to configure the generator toolbox" );
            return -1;
        }

        LMCINFO( lmclog, "Setting up simulation controller" );
        SimulationController controller;
        controller.SetFirstGenerator( toolbox.GetFirstGenerator() );
        if( ! controller.Configure( configurator.Config()->NodeAt( "simulation" ) ) )
        {
            LMCERROR( lmclog, "Unable to configure the simulation controller" );
            return -1;
        }

        LMCINFO( lmclog, "Preparing for run" );
        if( ! controller.Prepare() )
        {
            LMCERROR( lmclog, "Unable to prepare for the run" );
            return -1;
        }

        LMCINFO( lmclog, "Beginning simulation run" );
        controller.Run();

        LMCINFO( lmclog, "Run complete; finalizing" );
        controller.Finalize();
    }
    catch( Exception& e )
    {
        LMCERROR( lmclog, "Locust exception caught: " << e.what() );
        return -1;
    }
    catch( std::exception& e )
    {
        LMCERROR( lmclog, "Exception caught: " << e.what() );
        return -1;
    }

    return 0;
}

/*
 * LocustMC.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */


#include "application.hh"

#include "configurator.hh"
#include "LMCException.hh"
#include "LMCGeneratorToolbox.hh"
#include "LMCSimulationController.hh"

#include "logger.hh"

using namespace locust;
using namespace scarab;


LOGGER( lmclog, "LocustSim" );


void PrintHelpMessage()
{
    LPROG(lmclog, "\nUsage: LocustSim [options]\n\n" <<
           "  If using a config file, it should be specified as:  -c config_file.json\n" <<
           "  Config file options can be modified using:  \"address.of.option\"=value\n");
    return;
}


int main( int argc, char** argv )
{

	bool tScarabV1 = true;

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


        PrintHelpMessage();
        main_app the_main; // Scarab v2
        scarab::configurator configurator( argc, argv );  // Scarab v1

        LPROG( lmclog, "Setting up generator toolbox" );
        GeneratorToolbox toolbox;
        if( ! toolbox.Configure( configurator.config() ) )
        {
            LWARN( lmclog, "Unable to configure the v1 generator toolbox" );
            tScarabV1 = false;
        	CLI11_PARSE( the_main, argc, argv ); // Scarab v2
        	the_main.pre_callback(); // Scarab v2
        	if( ! toolbox.Configure( the_main.master_config().as_node() ))
        	{
        		LERROR( lmclog, "Unable to configure the v2 generator toolbox" );
        		return -1;
        	}
        }


        LPROG( lmclog, "Setting up simulation controller" );
        SimulationController controller;
        controller.SetFirstGenerator( toolbox.GetFirstGenerator() );

        if( tScarabV1 == true )
        {
        	if( ! controller.Configure( configurator.config()[ "simulation"].as_node()  ) )
        	{
                LERROR( lmclog, "Unable to configure the v1 simulation controller" );
                return -1;
        	}
        }
        else // Scarab v2
        {
            if( ! controller.Configure( the_main.master_config()["simulation"].as_node() ))
            	{
                LERROR( lmclog, "Unable to configure the v2 simulation controller" );
                return -1;
                }
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

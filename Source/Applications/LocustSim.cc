/*
 * LocustMC.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

//#include "application.hh"

#include "param.hh"
#include "param_json.hh"


//#include "configurator.hh"
#include "LMCException.hh"
#include "LMCGeneratorToolbox.hh"
#include "LMCSimulationController.hh"

#include "logger.hh"

using namespace locust;
using namespace scarab;

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

//        scarab::configurator configurator( argc, argv );
        //read file
        param_input_json t_input;
        param_ptr_t t_param_location( t_input.read_file( argv[1] ) );
        if( ! t_param_location )
        {
            LERROR( lmclog, "File did not read!" );
            return -1;
        }

        LINFO( lmclog, "File read and parsed: \n" << *t_param_location );

        //write file
        param_output_json t_output;
        bool t_did_write_file = t_output.write_file( *t_param_location, "test_output.json" );

        if( ! t_did_write_file )
        {
            LERROR( lmclog, "File did not write!" );
            return -1;
        }

        LINFO( lmclog, "File written successfully (test_output.json)" );

        LPROG( lmclog, "Setting up generator toolbox" );
        GeneratorToolbox toolbox;
//        if( ! toolbox.Configure( configurator.config() ) )
        if( ! toolbox.Configure( t_param_location->as_node() ) )
        {
            LERROR( lmclog, "Unable to configure the generator toolbox" );
            return -1;
        }


        LPROG( lmclog, "Setting up simulation controller" );
        SimulationController controller;
        controller.SetFirstGenerator( toolbox.GetFirstGenerator() );
        //        if( ! controller.Configure( configurator.config()[ "simulation"].as_node()  ) )

        if( ! controller.Configure( t_param_location->as_node() ))
        {
            LERROR( lmclog, "Unable to configure the simulation controller" );
            return -1;
        }

        LPROG( lmclog, "Successfully set up simulation controller" );

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

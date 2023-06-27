/*
 * LMCRunPause.cc
 *
 *  Created on: Jul 31, 2019
 *      Author: N.S. Oblath
 */

#include "LMCRunPause.hh"

#include "KSRun.h"
#include "logger.hh"



#include <csignal>

using namespace katrin;
namespace locust
{

    LOGGER( lmclog, "RunPause" );


    RunPause::RunPause() :
		    fToolbox(KToolbox::GetInstance()),
			KSComponent(),
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    RunPause::RunPause( const RunPause& aCopy ) :
		    fToolbox(KToolbox::GetInstance()),
            KSComponent(),
            fInterface( aCopy.fInterface )
    {
    }

    RunPause::~RunPause()
    {
    }

    RunPause* RunPause::Clone() const
    {
        return new RunPause( *this );
    }

    bool RunPause::ConfigureByInterface()
    {
    	if (fInterface->fConfigureKass)
    	{
    	    const scarab::param_node* aParam = fInterface->fConfigureKass->GetParameters();
    	    if (!this->Configure( *aParam ))
    	    {
    		    LERROR(lmclog,"Error configuring RunPause class");
    		    return false;
    	    }
    	}
    	else
    	{
		    LPROG(lmclog,"RunPause class did not need to be configured.");
		    return false;
    	}
        return true;
    }

    bool RunPause::Configure( const scarab::param_node& aParam )
     {

        if (fToolbox.IsInitialized())
        {
            for( auto sim : fToolbox.GetAll<Kassiopeia::KSSimulation>("project8_simulation"))
            {
            	// std::cout << "seed is " << sim->GetSeed();
                // TO-DO:  Change seed here, for pileup tests.
            }


            if( aParam.has( "waveguide-x" ) )
            {
            	// TO-DO:  Adjust Kass waveguide dimensions here.  Use fToolbox.

            }



            if( aParam.has( "cavity-radius" ) )
            {
            	if (fToolbox.HasKey("term_max_rlocust"))
            	{
                    fToolbox.Get<Kassiopeia::KSTermMaxR>("term_max_rlocust")->SetMaxR( aParam["cavity-radius"]().as_double() );
            	}
            	else
            	{
        		    LERROR(lmclog,"Can't find the terminator \"term_max_rlocust\" in the Kass xml file, which"
        		    		"means that electrons will not be properly terminated at the walls of the "
        		    		"Locust cavity with radius = \"cavity-radius\". ");
        		    exit(-1);
            	}
            }
        }
        else
        {
        	return false;
        }

    	return true;
     }



    bool RunPause::ExecutePreRunModification(Kassiopeia::KSRun &)
    {
    	ConfigureByInterface();

    	// Specifying this here because it comes before the interrupt in KSRoot.
    	fInterface->fKassEventReady = true;
        return true;
    }
//.
//.
// Between these two functions there will be an nevents interrupt in KSRoot.
//.
//.

    bool RunPause::ExecutePostRunModification(Kassiopeia::KSRun & aRun)
    {
    	//  No interrupt has happened yet in KSRoot.  Run still in progress.
//        fInterface->fRunInProgress = true;
        return true;
    }


} /* namespace locust */


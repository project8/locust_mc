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

        if ( fToolbox.IsInitialized() )
        {
            if( aParam.has( "electron-duration" ) )
             {
                 Kassiopeia::KSTermMaxTime* tLocustMaxTimeTerminator = new Kassiopeia::KSTermMaxTime();
                 tLocustMaxTimeTerminator->SetName("locust-electron-duration");
                 if (!fToolbox.HasKey("locust-electron-duration"))
                 {
                     tLocustMaxTimeTerminator->SetTime( aParam["electron-duration"]().as_double() );
                     tLocustMaxTimeTerminator->Initialize();
                     tLocustMaxTimeTerminator->Activate();
                     fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->AddTerminator(tLocustMaxTimeTerminator);
                 }
             }


            if( aParam.has( "cavity-radius" ) )
            {
                Kassiopeia::KSTermMaxR* tLocustMaxRTerminator = new Kassiopeia::KSTermMaxR();
                tLocustMaxRTerminator->SetName("locust-radius-terminator");
                if (!fToolbox.HasKey("locust-radius-terminator"))
                {
                    tLocustMaxRTerminator->SetMaxR( aParam["cavity-radius"]().as_double() );
                    tLocustMaxRTerminator->Initialize();
                    tLocustMaxRTerminator->Activate();
                    fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->AddTerminator(tLocustMaxRTerminator);
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
    	if ( !ConfigureByInterface() )
    	{
    	    return false;
    	}
    	else
    	{
    	    // Specifying this here because it comes before the interrupt in KSRoot.
    	    fInterface->fKassEventReady = true;
    	    return true;
    	}
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


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



            if (!fToolbox.HasKey("box_mockup_surface"))
            {

            auto* tBox = new KGeoBag::KGBoxSpace();
            tBox->XA(-4.e-3);
            tBox->XB(4.e-3);
            tBox->YA(-4.e-3);
            tBox->YB(4.e-3);
            tBox->ZA(-1.e-12);
            tBox->ZB(1.e-12);
            tBox->SetTag("box_mockup");

            auto* tSpace = new KGeoBag::KGSpace();
            tSpace->Volume(std::shared_ptr<KGeoBag::KGVolume>(tBox));

            fSurface = new Kassiopeia::KSGeoSurface();
            fSurface->SetName("box_mockup_surface");
            fSurface->AddContent(*tSpace->GetBoundaries()->begin());
            fToolbox.Add(fSurface);

            fLocustTermDeath = new Kassiopeia::KSTermDeath();
            fLocustTermDeath->Initialize(); // This is required.
//            fLocustTermDeath->Activate(); // not needed.
            fToolbox.Add(fLocustTermDeath);

            fCommand = fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->Command("add_terminator", fLocustTermDeath);
            fCommand->SetName("my_new_command");
//            tCommand->Activate(); // If this happens, the terminator is active all the time (bad).
            fToolbox.Get<Kassiopeia::KSGeoSpace>("space_world")->AddSurface(fSurface);
//            fSurface->Initialize();
//            fSurface->Activate();
            fSurface->AddCommand(fCommand);
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

/*
        fSurface->Deactivate();
        fSurface->Deinitialize();
    	fLocustTermDeath->Deinitialize();
    	fLocustTermDeath->Deactivate();
    	fSurface->RemoveCommand(fCommand);
    	fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->RemoveTerminator(fLocustTermDeath);
        fToolbox.Get<Kassiopeia::KSGeoSpace>("space_world")->RemoveSurface(fSurface);
*/

        return true;
    }


} /* namespace locust */


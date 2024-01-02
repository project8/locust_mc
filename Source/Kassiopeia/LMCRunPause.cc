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
                fLocustMaxTimeTerminator = new Kassiopeia::KSTermMaxTime();
                fLocustMaxTimeTerminator->SetName("locust-electron-duration");
                fLocustMaxTimeTerminator->SetTime( aParam["electron-duration"]().as_double() );
                fLocustMaxTimeTerminator->Initialize();
                fLocustMaxTimeTerminator->Activate();
                if (!fToolbox.HasKey("locust-electron-duration"))
                {
                    fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->AddTerminator(fLocustMaxTimeTerminator);
                }
            }

            if( aParam.has( "cavity-radius" ) )
            {
                fLocustMaxRTerminator = new Kassiopeia::KSTermMaxR();
                fLocustMaxRTerminator->SetName("locust-radius-terminator");
                fLocustMaxRTerminator->SetMaxR( aParam["cavity-radius"]().as_double() );
                fLocustMaxRTerminator->Initialize();
                fLocustMaxRTerminator->Activate();
                if (!fToolbox.HasKey("locust-radius-terminator"))
                {
                    fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->AddTerminator(fLocustMaxRTerminator);
                }
            }

            if ( aParam.has( "waveguide-x" ) )
            {
                if ( aParam.has( "waveguide-y" ) && aParam.has( "waveguide-z" ) )
                {
                    fBox = new KGeoBag::KGBoxSpace();
                    fBox->XA(-aParam["waveguide-x"]().as_double()/2.);
                    fBox->XB(aParam["waveguide-x"]().as_double()/2.);
                    fBox->YA(-aParam["waveguide-y"]().as_double()/2.);
                    fBox->YB(aParam["waveguide-y"]().as_double()/2.);
                    fBox->ZA(-aParam["waveguide-z"]().as_double()/2.);
                    fBox->ZB(aParam["waveguide-z"]().as_double()/2.);
                    fBox->SetTag("waveguide_box");

                    fKGSpace = new KGeoBag::KGSpace();
                    fKGSpace->Volume(std::shared_ptr<KGeoBag::KGVolume>(fBox));

                    fSurface = new Kassiopeia::KSGeoSurface();
                    fSurface->SetName("waveguide_surfaces");
                    auto it = fKGSpace->GetBoundaries()->begin();
                    while (it != fKGSpace->GetBoundaries()->end())
                    {
                        fSurface->AddContent(*it);
                        *it++;
                    }

                    if (!fToolbox.HasKey("waveguide_surfaces"))
                    {
                        fToolbox.Add(fSurface);

                        fLocustTermDeath = new Kassiopeia::KSTermDeath();
                        fLocustTermDeath->Initialize();
                        fToolbox.Add(fLocustTermDeath);

                        fCommand = fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->Command("add_terminator", fLocustTermDeath);
                        fKSSpaces = fToolbox.GetAll<Kassiopeia::KSGeoSpace>();
                        if ( fKSSpaces.size() == 1 )
                        {
                            LPROG(lmclog,"LMCRunPause found the KSGeoSpace named <" << fKSSpaces[0]->GetName() << ">");
                            fKSSpaces[0]->AddSurface(fSurface);
                        }
                        else
                        {
                            LERROR(lmclog,"Only one KSGeoSpace instance was expected to be in the KToolbox.");
                            exit(-1);
                        }

                        fSurface->AddCommand(fCommand);
                    }
                }
                else
                {
        		    LERROR(lmclog,"Waveguide dimensions waveguide-x, waveguide-y, and waveguide-z are needed.");
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

    bool RunPause::DeleteLocalKassObjects()
    {

        delete fLocustMaxTimeTerminator;
        delete fLocustMaxRTerminator;
        delete fBox;
        delete fKGSpace;
        delete fSurface;
        delete fLocustTermDeath;
        delete fCommand;
        for (unsigned i=0; i<fKSSpaces.size(); i++)
        {
        	delete fKSSpaces[0];
        }


    	return true;
    }

    bool RunPause::ExecutePostRunModification(Kassiopeia::KSRun & aRun)
    {
    	//  No interrupt has happened yet in KSRoot.  Run still in progress.
//        fInterface->fRunInProgress = true;

        if ( !fToolbox.IsInitialized() )
        {
        	DeleteLocalKassObjects();
        }

        return true;
    }


} /* namespace locust */


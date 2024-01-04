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


/*
            auto tGen = fToolbox.GetAll<Kassiopeia::KSGenerator>();
            for (unsigned i=0; i<tGen.size(); i++)
            {
            	if ( tGen[i]->IsActivated() )
            	{
                    fToolbox.Get<Kassiopeia::KSRootGenerator>("root_generator")->ClearGenerator(tGen[i]);
            	}
            }

            Kassiopeia::KSGenEnergyComposite* tGenEnergyComposite;
            Kassiopeia::KSGenValueUniform* tEnergyGenerator = new Kassiopeia::KSGenValueUniform();
            tEnergyGenerator->SetValueMin(18600.);
            tEnergyGenerator->SetValueMax(18600.);
            tGenEnergyComposite->SetEnergyValue(tEnergyGenerator);

            Kassiopeia::KSGenPositionRectangularComposite* tGenPositionComposite = new Kassiopeia::KSGenPositionRectangularComposite();
            tGenPositionComposite->SetOrigin(GetKSWorldSpace()->GetContent()->)
            Kassiopeia::KSGenValueUniform* tPositionXGenerator = new Kassiopeia::KSGenValueUniform();
            Kassiopeia::KSGenValueUniform* tPositionYGenerator = new Kassiopeia::KSGenValueUniform();
            Kassiopeia::KSGenValueUniform* tPositionZGenerator = new Kassiopeia::KSGenValueUniform();
            tPositionXGenerator->SetValueMin(0.004);
            tPositionXGenerator->SetValueMin(0.004);
            tPositionYGenerator->SetValueMin(0.0);
            tPositionYGenerator->SetValueMin(0.0);
            tPositionZGenerator->SetValueMin(0.0);
            tPositionZGenerator->SetValueMin(0.0);
            tGenPositionComposite->SetXValue(tPositionXGenerator);
            tGenPositionComposite->SetYValue(tPositionYGenerator);
            tGenPositionComposite->SetZValue(tPositionZGenerator);


            tPositionGenerator->SetXValue()
//            fGenerator->

*/


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

                    fKGSpace = GetKGWorldSpace();
                    KGeoBag::KGSpace* tKGSpace = new KGeoBag::KGSpace();
                    tKGSpace->Volume(std::shared_ptr<KGeoBag::KGVolume>(fBox));
                    fKGSpace->GetChildSpaces()->at(0)->AddChildSpace(tKGSpace);

                    fSurface = new Kassiopeia::KSGeoSurface();
                    fSurface->SetName("waveguide_surfaces");
                    auto it = tKGSpace->GetBoundaries()->begin();
                    while (it != tKGSpace->GetBoundaries()->end())
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
                        fKSSpace = GetKSWorldSpace();
                        fKSSpace->AddSurface(fSurface);
                        fSurface->AddCommand(fCommand);

                        KGeoBag::KGSpace* test = GetKGWorldSpace();
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



    Kassiopeia::KSGeoSpace* RunPause::GetKSWorldSpace()
    {
    	std::vector<Kassiopeia::KSGeoSpace*> tKSSpaces = fToolbox.GetAll<Kassiopeia::KSGeoSpace>();
        if ( tKSSpaces.size() == 1 )
        {
            LPROG(lmclog,"LMCRunPause found the KSGeoSpace named <" << tKSSpaces[0]->GetName() << ">");
        }
        else
        {
            LERROR(lmclog,"One and only one KSGeoSpace instance was expected to be in the KToolbox.");
            getchar();
            exit(-1);
        }

        return tKSSpaces[0];
    }

    KGeoBag::KGSpace* RunPause::GetKGWorldSpace()
    {
    	Kassiopeia::KSGeoSpace* tKSSpace = GetKSWorldSpace();
    	std::vector<KGeoBag::KGSpace*> tKGSpaces = tKSSpace->GetContent();
        if (( tKGSpaces.size() == 1 ) && (tKGSpaces[0]->GetChildSpaces()->at(0)->GetName()=="project8"))
        {
            LPROG(lmclog,"LMCRunPause found the KGSpace named <" << tKGSpaces[0]->GetName() << "/project8>");
        }
        else
        {
            LERROR(lmclog,"One and only one KGSpace instance was expected to be in the KToolbox, and "
            		"it was expected to contain a child KGSpace named <project8>");
            getchar();
            exit(-1);
        }

        const KGeoBag::KGSpace* test = tKGSpaces[0];

        return tKGSpaces[0];
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
        delete fKSSpace;

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


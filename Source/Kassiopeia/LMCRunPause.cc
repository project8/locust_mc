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
            return true;
    	}
        return true;
    }

    bool RunPause::Configure( const scarab::param_node& aParam )
     {

        if ( fToolbox.IsInitialized() )
        {
            if( aParam.has( "track-length" ) || ( aParam.has( "random-track-length" ) ) )
            {
                if (!AddMaxTimeTerminator( aParam ))
                {
                    return false;
                }
            }

            if( aParam.has( "cavity-radius" ) )
            {
                if (!AddMaxRTerminator( aParam ))
                {
                    return false;
                }
            }

            if ( aParam.has( "waveguide-x" ) )
            {
                if (!AddWaveguideTerminator( aParam ))
                {
                    return false;
                }
            }

            if ( aParam.has( "ks-starting-xpos-min" )  )
            {
                if (!AddGenerator( aParam ))
                {
                    return false;
                }
            }

            return true;

        }
        else
        {
            return false;
        }

     }

    bool RunPause::AddGenerator( const scarab::param_node& aParam )
    {

        if ( aParam.has( "ks-starting-xpos-min" ) && aParam.has( "ks-starting-xpos-max" )
        &&   aParam.has( "ks-starting-ypos-min" ) && aParam.has( "ks-starting-ypos-max" )
		&&   aParam.has( "ks-starting-zpos-min" ) && aParam.has( "ks-starting-zpos-max" )
		&&   aParam.has( "ks-starting-energy-min" ) && aParam.has( "ks-starting-energy-max" )
		&&   aParam.has( "ks-starting-pitch-min" ) && aParam.has( "ks-starting-pitch-max" ) )
        {
            KRandom::GetInstance().SetSeed(time(NULL));
            auto tGen = fToolbox.GetAll<Kassiopeia::KSGenerator>();
            for (unsigned i=0; i<tGen.size(); i++)
            {
                if ( (tGen[i]->IsActivated()) && (tGen[i]->GetName()!="root_generator") )
                {
                    LPROG(lmclog,"Clearing " << tGen[i]->GetName());
                    fToolbox.Get<Kassiopeia::KSRootGenerator>("root_generator")->ClearGenerator(tGen[i]);
                    fToolbox.Remove(tGen[i]->GetName());
                }
            }

            Kassiopeia::KSGenValueFix* tGenPidComposite = new Kassiopeia::KSGenValueFix();
            tGenPidComposite->SetValue(11); // electron

            Kassiopeia::KSGenTimeComposite* tGenTimeComposite = new Kassiopeia::KSGenTimeComposite();
            Kassiopeia::KSGenValueUniform* tTimeGenerator = new Kassiopeia::KSGenValueUniform();
            tTimeGenerator->SetValueMin(0.);
            tTimeGenerator->SetValueMax(0.);
            tGenTimeComposite->SetTimeValue(tTimeGenerator);

            Kassiopeia::KSGenEnergyComposite* tGenEnergyComposite = new Kassiopeia::KSGenEnergyComposite();
            Kassiopeia::KSGenValueUniform* tEnergyGenerator = new Kassiopeia::KSGenValueUniform();
            if ( aParam.has( "ks-starting-energy-min" ) && ( aParam.has( "ks-starting-energy-max" ) ) )
            {
                tEnergyGenerator->SetValueMin( aParam["ks-starting-energy-min"]().as_double() ); // eV
                tEnergyGenerator->SetValueMax( aParam["ks-starting-energy-max"]().as_double() ); // eV
            }
            else
            {
                tEnergyGenerator->SetValueMin( 18600. ); // eV
                tEnergyGenerator->SetValueMax( 18600. ); // eV
            }
            tGenEnergyComposite->SetEnergyValue(tEnergyGenerator);

            Kassiopeia::KSGenPositionRectangularComposite* tGenPositionComposite = new Kassiopeia::KSGenPositionRectangularComposite();
            tGenPositionComposite->SetOrigin(GetKGWorldSpace()->GetOrigin());
            Kassiopeia::KSGenValueUniform* tPositionXGenerator = new Kassiopeia::KSGenValueUniform();
            Kassiopeia::KSGenValueUniform* tPositionYGenerator = new Kassiopeia::KSGenValueUniform();
            Kassiopeia::KSGenValueUniform* tPositionZGenerator = new Kassiopeia::KSGenValueUniform();
            tPositionXGenerator->SetValueMin( aParam["ks-starting-xpos-min"]().as_double() ); // meters
            tPositionXGenerator->SetValueMax( aParam["ks-starting-xpos-max"]().as_double() );
            tPositionYGenerator->SetValueMin( aParam["ks-starting-ypos-min"]().as_double() );
            tPositionYGenerator->SetValueMax( aParam["ks-starting-ypos-max"]().as_double() );
            tPositionZGenerator->SetValueMin( aParam["ks-starting-zpos-min"]().as_double() );
            tPositionZGenerator->SetValueMax( aParam["ks-starting-zpos-max"]().as_double() );
            tGenPositionComposite->SetXValue(tPositionXGenerator);
            tGenPositionComposite->SetYValue(tPositionYGenerator);
            tGenPositionComposite->SetZValue(tPositionZGenerator);

            Kassiopeia::KSGenDirectionSphericalComposite* tGenDirectionComposite = new Kassiopeia::KSGenDirectionSphericalComposite();
            Kassiopeia::KSGenValueUniform* tThetaGenerator = new Kassiopeia::KSGenValueUniform();
            tThetaGenerator->SetValueMin( aParam["ks-starting-pitch-min"]().as_double() );
            tThetaGenerator->SetValueMax( aParam["ks-starting-pitch-max"]().as_double() );
            Kassiopeia::KSGenValueUniform* tPhiGenerator = new Kassiopeia::KSGenValueUniform();
            if ( aParam.has( "ks-starting-phi-min" ) && ( aParam.has( "ks-starting-phi-max" ) ) )
            {
                tPhiGenerator->SetValueMin( aParam["ks-starting-phi-min"]().as_double() );
                tPhiGenerator->SetValueMax( aParam["ks-starting-phi-max"]().as_double() );
            }
            else
            {
                tPhiGenerator->SetValueMin( 0. );
                tPhiGenerator->SetValueMax( 0. );
            }
            tGenDirectionComposite->SetPhiValue(tPhiGenerator);
            tGenDirectionComposite->SetThetaValue(tThetaGenerator);

            fGenerator = new Kassiopeia::KSGenGeneratorComposite();
            fGenerator->SetPid(tGenPidComposite);
            fGenerator->AddCreator(tGenPositionComposite);
            fGenerator->AddCreator(tGenEnergyComposite);
            fGenerator->AddCreator(tGenDirectionComposite);
            fGenerator->AddCreator(tGenTimeComposite);
            fGenerator->SetName("gen_project8");
            fGenerator->Initialize();
            fGenerator->Activate();

            if (!fToolbox.HasKey("gen_project8"))
            {
                fToolbox.Add(fGenerator);
                fToolbox.Get<Kassiopeia::KSRootGenerator>("root_generator")->SetGenerator(fGenerator);
            }

        }

        else
        {
            LERROR(lmclog,"To configure starting e- kinematics, all of these parameters are needed:  "
            	"ks-starting-xpos-min, ks-starting-xpos-max, ks-starting-ypos-min, ks-starting-ypos-max, "
            	"ks-starting-zpos-min, ks-starting-zpos-max, ks-starting-pitch-min, ks-starting-pitch-max,"
            	"ks-starting-energy-min, ks-starting-energy-max ");
            return false;
        }

    	return true;
    }


    bool RunPause::AddMaxRTerminator( const scarab::param_node& aParam )
    {
    	/* Remove any existing KSTermMaxR objects */
        auto tMaxR = fToolbox.GetAll<Kassiopeia::KSTermMaxR>();
        for (unsigned i=0; i<tMaxR.size(); i++)
        {
        	fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->RemoveTerminator(tMaxR[i]);
        	fToolbox.Remove(tMaxR[i]->GetName());
        }

        if (!fToolbox.HasKey("ksmax-r-project8"))
        {
            fLocustMaxRTerminator = new Kassiopeia::KSTermMaxR();
            fLocustMaxRTerminator->SetName("ksmax-r-project8");
            fLocustMaxRTerminator->SetMaxR( aParam["cavity-radius"]().as_double() );
            fLocustMaxRTerminator->Initialize();
            fToolbox.Add(fLocustMaxRTerminator);
            fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->AddTerminator(fLocustMaxRTerminator);
        }
        return true;
    }


    int RunPause::GetSeed( const scarab::param_node& aParam )
    {
        int tSeed = 0;
        if ( aParam.has( "random-track-seed" ) )
        {
            tSeed = aParam["random-track-seed"]().as_int();
        }
        else
        {
            tSeed = time(NULL);
        }
        LPROG(lmclog,"Setting random seed for track length to " << tSeed);
        return tSeed;
    }


    bool RunPause::AddMaxTimeTerminator( const scarab::param_node& aParam )
    {

    	/* Remove any existing KSTermMaxTime objects */
        auto tMaxTime = fToolbox.GetAll<Kassiopeia::KSTermMaxTime>();
        for (unsigned i=0; i<tMaxTime.size(); i++)
        {
        	fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->RemoveTerminator(tMaxTime[i]);
        	fToolbox.Remove(tMaxTime[i]->GetName());
        }

        if (!fToolbox.HasKey("ksmax-time-project8"))
        {

            double tMaxTrackLength = 0.;
            fLocustMaxTimeTerminator = new Kassiopeia::KSTermMaxTime();

    	    if ( aParam.has( "track-length" ) )
    	    {
    	        tMaxTrackLength = aParam["track-length"]().as_double();
    	        fLocustMaxTimeTerminator->SetTime( tMaxTrackLength );
    	    }
    	    else
    	    {
                LERROR(lmclog,"Parameter \"track-length\" is required to set the maximum track length.");
                return false;
    	    }

    	    if ( aParam.has( "random-track-length" ) )
    	    {
                if ( aParam["random-track-length"]().as_bool() == true)
                {
                    srand ( GetSeed( aParam ));
                    double tRandomTime = tMaxTrackLength/10. * ( 1 + rand() % 10 ); // 0.1*tMaxTrackLength < t < 1.1*tMaxTrackLength
                    fLocustMaxTimeTerminator->SetTime( tRandomTime );
                    LPROG(lmclog,"Randomizing the track length to " << tRandomTime);
                }
    	    }

            fLocustMaxTimeTerminator->SetName("ksmax-time-project8");
            fLocustMaxTimeTerminator->Initialize();
            fLocustMaxTimeTerminator->Activate();
            fToolbox.Add(fLocustMaxTimeTerminator);
            fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->AddTerminator(fLocustMaxTimeTerminator);
        }
        return true;
    }

    bool RunPause::AddWaveguideTerminator( const scarab::param_node& aParam )
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
            fSurface->SetName("waveguide_surfaces_project8");
            auto it = tKGSpace->GetBoundaries()->begin();
            while (it != tKGSpace->GetBoundaries()->end())
            {
                fSurface->AddContent(*it);
                *it++;
            }

            if (!fToolbox.HasKey("waveguide_surfaces_project8"))
            {
                fToolbox.Add(fSurface);

                fLocustTermDeath = new Kassiopeia::KSTermDeath();
                fLocustTermDeath->Initialize();
                fToolbox.Add(fLocustTermDeath);

                fCommand = fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->Command("add_terminator", fLocustTermDeath);
                fKSSpace = GetKSWorldSpace();
                fKSSpace->AddSurface(fSurface);
                fSurface->AddCommand(fCommand);
            }
            return true;
        }
        else
        {
            LERROR(lmclog,"Waveguide dimensions waveguide-x, waveguide-y, and waveguide-z are needed.");
            return false;
        }

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
            LERROR(lmclog,"One and only one KGSpace instance was expected to be in the KToolbox, and/or "
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

        fLocustMaxTimeTerminator = NULL;
        fLocustMaxRTerminator = NULL;
        fBox = NULL;
        fKGSpace = NULL;
        fSurface = NULL;
        fLocustTermDeath = NULL;
        fCommand = NULL;
        fKSSpace = NULL;
        fGenerator = NULL;

    	return true;
    }

    bool RunPause::ExecutePostRunModification(Kassiopeia::KSRun & aRun)
    {
    	//  No interrupt has happened yet in KSRoot.  Run still in progress.
//        fInterface->fRunInProgress = true;


    	DeleteLocalKassObjects();


        return true;
    }


} /* namespace locust */


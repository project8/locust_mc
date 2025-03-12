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
        fLocustMaxTimeTerminator( nullptr ),
        fLocustMaxRTerminator( nullptr ),
        fBox( nullptr ),
        fKGSpace( nullptr ),
        fSurface( nullptr ),
        fLocustTermDeath( nullptr ),
        fCommand( nullptr ),
        fKSSpace( nullptr ),
		fGenDirectionComposite( nullptr ),
		fThetaGenerator( nullptr ),
		fPhiGenerator( nullptr ),
		fGenPositionComposite( nullptr ),
		fPositionXGenerator( nullptr ),
		fPositionYGenerator( nullptr ),
		fPositionZGenerator( nullptr ),
		fGenEnergyComposite( nullptr ),
		fEnergyGenerator( nullptr ),
		fTimeGenerator( nullptr ),
		fGenTimeComposite( nullptr ),
		fGenPidComposite( nullptr ),
        fGenerator( nullptr ),
        fMinTrackLengthFraction(0.1),
        fConfigurationComplete( false ),
        fEventCounter( 0 ),
        fMaxEvents( 1 ),
        fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    RunPause::RunPause( const RunPause& aCopy ) :
        fToolbox(KToolbox::GetInstance()),
        KSComponent(),
        fLocustMaxTimeTerminator( nullptr ),
        fLocustMaxRTerminator( nullptr ),
        fBox( nullptr ),
        fKGSpace( nullptr ),
        fSurface( nullptr ),
        fLocustTermDeath( nullptr ),
        fCommand( nullptr ),
        fKSSpace( nullptr ),
		fGenDirectionComposite( nullptr ),
		fThetaGenerator( nullptr ),
		fPhiGenerator( nullptr ),
		fGenPositionComposite( nullptr ),
		fPositionXGenerator( nullptr ),
		fPositionYGenerator( nullptr ),
		fPositionZGenerator( nullptr ),
		fGenEnergyComposite( nullptr ),
		fEnergyGenerator( nullptr ),
		fTimeGenerator( nullptr ),
		fGenTimeComposite( nullptr ),
		fGenPidComposite( nullptr),
		fGenerator( nullptr ),
        fMinTrackLengthFraction(0.1),
        fConfigurationComplete( false ),
        fEventCounter( 0 ),
        fMaxEvents( 1 ),
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
        if ( !fConfigurationComplete )
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

            fConfigurationComplete = true;
            LPROG(lmclog,"RunPause has been configured.");
        }

        return true;
    }

    bool RunPause::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "n-max-events" ) )
        {
            fMaxEvents = aParam["n-max-events"]().as_double();
        }

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

        if (!fToolbox.HasKey("gen_project8"))
        {

            if ( aParam.has( "ks-starting-xpos-min" ) && aParam.has( "ks-starting-xpos-max" )
            &&   aParam.has( "ks-starting-ypos-min" ) && aParam.has( "ks-starting-ypos-max" )
		    &&   aParam.has( "ks-starting-zpos-min" ) && aParam.has( "ks-starting-zpos-max" )
		    &&   aParam.has( "ks-starting-energy-min" ) && aParam.has( "ks-starting-energy-max" )
		    &&   aParam.has( "ks-starting-pitch-min" ) && aParam.has( "ks-starting-pitch-max" ) )
            {
                int tSeed = 0;
                if ( aParam.has( "random-track-seed" ) )
                {
                    tSeed = aParam["random-track-seed"]().as_int();
                	KRandom::GetInstance().SetSeed(tSeed);
                }
                else
                {
                    GetSeed( aParam );
                }

                auto tGen = fToolbox.GetAll<Kassiopeia::KSGenerator>();
                for (unsigned i=0; i<tGen.size(); i++)
                {
                    if ( (tGen[i]->IsActivated()) &&  (tGen[i]->GetName()!="root_generator") )
                    {
                        LPROG(lmclog,"Clearing " << tGen[i]->GetName() << " from KSRoot ... ");
                        fToolbox.Get<Kassiopeia::KSRootGenerator>("root_generator")->ClearGenerator(tGen[i]);
                        tGen[i]->RemoveTags(tGen[i]->GetTags());
                        tGen[i]->Deactivate();
                        tGen[i]->Deinitialize();
                        fToolbox.Remove(tGen[i]->GetName());
                    }
                }

                if ( fGenPidComposite == nullptr ) fGenPidComposite = new Kassiopeia::KSGenValueFix();
                fGenPidComposite->SetValue(11); // electron

                if ( fGenTimeComposite == nullptr ) fGenTimeComposite = new Kassiopeia::KSGenTimeComposite();
                if ( fTimeGenerator == nullptr ) fTimeGenerator = new Kassiopeia::KSGenValueUniform();
                fTimeGenerator->SetValueMin(0.);
                fTimeGenerator->SetValueMax(0.);
                fGenTimeComposite->SetTimeValue(fTimeGenerator);

                if ( fGenEnergyComposite == nullptr ) fGenEnergyComposite = new Kassiopeia::KSGenEnergyComposite();
                if ( fEnergyGenerator == nullptr ) fEnergyGenerator = new Kassiopeia::KSGenValueUniform();
                if ( aParam.has( "ks-starting-energy-min" ) && ( aParam.has( "ks-starting-energy-max" ) ) )
                {
                    fEnergyGenerator->SetValueMin( aParam["ks-starting-energy-min"]().as_double() ); // eV
                    fEnergyGenerator->SetValueMax( aParam["ks-starting-energy-max"]().as_double() ); // eV
                }
                else
                {
                    fEnergyGenerator->SetValueMin( 18600. ); // eV
                    fEnergyGenerator->SetValueMax( 18600. ); // eV
                }

                fGenEnergyComposite->SetEnergyValue(fEnergyGenerator);

                if ( fGenPositionComposite == nullptr ) fGenPositionComposite = new Kassiopeia::KSGenPositionRectangularComposite();
                fGenPositionComposite->SetOrigin(GetKGWorldSpace()->GetOrigin());
                if ( fPositionXGenerator == nullptr ) fPositionXGenerator = new Kassiopeia::KSGenValueUniform();
                if ( fPositionYGenerator == nullptr ) fPositionYGenerator = new Kassiopeia::KSGenValueUniform();
                if ( fPositionZGenerator == nullptr ) fPositionZGenerator = new Kassiopeia::KSGenValueUniform();
                fPositionXGenerator->SetValueMin( aParam["ks-starting-xpos-min"]().as_double() ); // meters
                fPositionXGenerator->SetValueMax( aParam["ks-starting-xpos-max"]().as_double() );
                fPositionYGenerator->SetValueMin( aParam["ks-starting-ypos-min"]().as_double() );
                fPositionYGenerator->SetValueMax( aParam["ks-starting-ypos-max"]().as_double() );
                fPositionZGenerator->SetValueMin( aParam["ks-starting-zpos-min"]().as_double() );
                fPositionZGenerator->SetValueMax( aParam["ks-starting-zpos-max"]().as_double() );
                fGenPositionComposite->SetXValue(fPositionXGenerator);
                fGenPositionComposite->SetYValue(fPositionYGenerator);
                fGenPositionComposite->SetZValue(fPositionZGenerator);

                if ( fGenDirectionComposite == nullptr ) fGenDirectionComposite = new Kassiopeia::KSGenDirectionSphericalComposite();
                if ( fThetaGenerator == nullptr ) fThetaGenerator = new Kassiopeia::KSGenValueAngleSpherical();
                fThetaGenerator->SetAngleMin( aParam["ks-starting-pitch-min"]().as_double() );
                fThetaGenerator->SetAngleMax( aParam["ks-starting-pitch-max"]().as_double() );
                if ( fPhiGenerator == nullptr ) fPhiGenerator = new Kassiopeia::KSGenValueUniform();
                if ( aParam.has( "ks-starting-phi-min" ) && ( aParam.has( "ks-starting-phi-max" ) ) )
                {
                    fPhiGenerator->SetValueMin( aParam["ks-starting-phi-min"]().as_double() );
                    fPhiGenerator->SetValueMax( aParam["ks-starting-phi-max"]().as_double() );
                }
                else
                {
                    fPhiGenerator->SetValueMin( 0. );
                    fPhiGenerator->SetValueMax( 0. );
                }
                fGenDirectionComposite->SetPhiValue(fPhiGenerator);
                fGenDirectionComposite->SetThetaValue(fThetaGenerator);

                if ( fGenerator == nullptr ) fGenerator = new Kassiopeia::KSGenGeneratorComposite();

                fGenerator->SetPid(fGenPidComposite);
                fGenerator->AddCreator(fGenPositionComposite);
                fGenerator->AddCreator(fGenEnergyComposite);
                fGenerator->AddCreator(fGenDirectionComposite);
                fGenerator->AddCreator(fGenTimeComposite);
                fGenerator->SetName("gen_project8");
                fGenerator->Initialize();
                fGenerator->Activate();

                fToolbox.Add(fGenerator);
                fToolbox.Get<Kassiopeia::KSRootGenerator>("root_generator")->SetGenerator(fGenerator);
                LPROG(lmclog,"\"gen-project8\" has just been added to the KToolbox.");
            }

            else
            {
                LERROR(lmclog,"To configure starting e- kinematics, all of these parameters are needed:  "
            	    "ks-starting-xpos-min, ks-starting-xpos-max, ks-starting-ypos-min, ks-starting-ypos-max, "
            	    "ks-starting-zpos-min, ks-starting-zpos-max, ks-starting-pitch-min, ks-starting-pitch-max,"
            	    "ks-starting-energy-min, ks-starting-energy-max ");
                return false;
            }
        }
        else
        {
            LPROG(lmclog,"\"gen-project8\" is already in the KToolbox.");
        }

    	return true;
    }


    bool RunPause::AddMaxRTerminator( const scarab::param_node& aParam )
    {
        if (!fToolbox.HasKey("ksmax-r-project8"))
        {
            /* Remove any existing KSTermMaxR objects */
            auto tMaxR = fToolbox.GetAll<Kassiopeia::KSTermMaxR>();
            for (unsigned i=0; i<tMaxR.size(); i++)
            {
        	    fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->RemoveTerminator(tMaxR[i]);
        	    fToolbox.Remove(tMaxR[i]->GetName());
            }

            if ( fLocustMaxRTerminator == nullptr ) fLocustMaxRTerminator = new Kassiopeia::KSTermMaxR();
            fLocustMaxRTerminator->SetName("ksmax-r-project8");
            fLocustMaxRTerminator->SetMaxR( aParam["cavity-radius"]().as_double() );
            fLocustMaxRTerminator->Initialize();
            fLocustMaxRTerminator->Activate();
            fToolbox.Add(fLocustMaxRTerminator);
            fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->AddTerminator(fLocustMaxRTerminator);
            LPROG(lmclog,"\"ksmax-r-project8\" has just been added to the KToolbox.");
        }
        else
        {
            LPROG(lmclog,"\"ksmax-r-project8\" is already in the KToolbox.");
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
            struct timeval tv;
            gettimeofday(&tv, NULL);
            tSeed = tv.tv_usec;
        }
        LPROG(lmclog,"Setting random seed for track length to " << tSeed);
        return tSeed;
    }


    bool RunPause::AddMaxTimeTerminator( const scarab::param_node& aParam )
    {
        if (!fToolbox.HasKey("ksmax-time-project8"))
        {

            /* Remove any existing KSTermMaxTime objects */
            auto tMaxTime = fToolbox.GetAll<Kassiopeia::KSTermMaxTime>();
            for (unsigned i=0; i<tMaxTime.size(); i++)
            {
                fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->RemoveTerminator(tMaxTime[i]);
                fToolbox.Remove(tMaxTime[i]->GetName());
            }

            double tMaxTrackLength = 0.;
            if ( fLocustMaxTimeTerminator == nullptr ) fLocustMaxTimeTerminator = new Kassiopeia::KSTermMaxTime();

    	    if ( aParam.has( "min-track-length-fraction" ) )
            {
                fMinTrackLengthFraction = aParam["min-track-length-fraction"]().as_double();
                LPROG(lmclog,"Setting minimum track length fraction to " << fMinTrackLengthFraction);
    	    }

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
                    scarab::param_node default_setting;
                    default_setting.add("name","uniform");
                    fTrackLengthDistribution = fDistributionInterface.get_dist(default_setting);
                    fDistributionInterface.SetSeed( GetSeed(aParam) );
                    double tMinTrackLength = tMaxTrackLength * fMinTrackLengthFraction;
                    double tRandomTime = tMinTrackLength + (tMaxTrackLength - tMinTrackLength) * fTrackLengthDistribution->Generate();
                    fLocustMaxTimeTerminator->SetTime( tRandomTime );
                    LPROG(lmclog,"Randomizing the track length to " << tRandomTime);
                }
            }

            fLocustMaxTimeTerminator->SetName("ksmax-time-project8");
            fLocustMaxTimeTerminator->Initialize();
            fLocustMaxTimeTerminator->Activate();
            fToolbox.Add(fLocustMaxTimeTerminator);
            fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->AddTerminator(fLocustMaxTimeTerminator);
            LPROG(lmclog,"\"ksmax-time-project8\" has just been added to the KToolbox.");
        }
        else
        {
            LPROG(lmclog,"\"ksmax-time-project8\" is already in the KToolbox.");
        }
        return true;
    }

    bool RunPause::AddWaveguideTerminator( const scarab::param_node& aParam )
    {
        if (!fToolbox.HasKey("waveguide_surfaces_project8"))
        {

            if ( aParam.has( "waveguide-y" ) && aParam.has( "waveguide-z" ) )
            {
                if ( fBox == nullptr ) fBox = new KGeoBag::KGBoxSpace();
                fBox->XA(-aParam["waveguide-x"]().as_double()/2.);
                fBox->XB(aParam["waveguide-x"]().as_double()/2.);
                fBox->YA(-aParam["waveguide-y"]().as_double()/2.);
                fBox->YB(aParam["waveguide-y"]().as_double()/2.);
                fBox->ZA(-aParam["waveguide-z"]().as_double()/2.);
                fBox->ZB(aParam["waveguide-z"]().as_double()/2.);
                fBox->SetTag("waveguide_box");

                if ( fKGSpace == nullptr ) fKGSpace = GetKGWorldSpace();
                KGeoBag::KGSpace* tKGSpace = new KGeoBag::KGSpace();
                tKGSpace->Volume(std::shared_ptr<KGeoBag::KGVolume>(fBox));
                fKGSpace->GetChildSpaces()->at(0)->AddChildSpace(tKGSpace);

                if ( fSurface == nullptr ) fSurface = new Kassiopeia::KSGeoSurface();
                fSurface->SetName("waveguide_surfaces_project8");
                auto it = tKGSpace->GetBoundaries()->begin();
                while (it != tKGSpace->GetBoundaries()->end())
                {
                    fSurface->AddContent(*it);
                    *it++;
                }

                fToolbox.Add(fSurface);
                fLocustTermDeath = new Kassiopeia::KSTermDeath();
                fLocustTermDeath->Initialize();
                fLocustTermDeath->Activate();
                fToolbox.Add(fLocustTermDeath);
                fLocustTermDeath->SetName("waveguide_death");
                fCommand = fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->Command("add_terminator", fLocustTermDeath);
                if ( fKSSpace == nullptr ) fKSSpace = GetKSWorldSpace();
                fKSSpace->AddSurface(fSurface);
                fSurface->AddCommand(fCommand);
            }
        }
        else
        {
            LERROR(lmclog,"Waveguide dimensions waveguide-x, waveguide-y, and waveguide-z are needed.");
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

    bool RunPause::DeleteLocalKassObjects()
    {
        // fMaxTimeTerminator
        if (fToolbox.HasKey("ksmax-time-project8"))
        {
            fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->RemoveTerminator(fLocustMaxTimeTerminator);
            fLocustMaxTimeTerminator = NULL;
            LPROG(lmclog,"Removing fLocustMaxTimeTerminator from KSRoot ...");
        }

        // fMaxRTerminator
        if (fToolbox.HasKey("ksmax-r-project8"))
        {
            fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->RemoveTerminator(fLocustMaxRTerminator);
            fLocustMaxRTerminator = NULL;
            LPROG(lmclog,"Removing fLocustMaxRTerminator from KSRoot ...");
        }

        // fSurface
        if (fToolbox.HasKey("waveguide_surfaces_project8"))
        {
            fToolbox.Get<Kassiopeia::KSRootTerminator>("root_terminator")->RemoveTerminator(fLocustTermDeath);
            fLocustTermDeath = NULL;
            fCommand = NULL;
            fBox = NULL;
            fSurface = NULL;
            LPROG(lmclog,"Removing fSurface from KSRoot ...");
        }

        // fGenerator
        if (fToolbox.HasKey("gen_project8"))
        {
            fToolbox.Get<Kassiopeia::KSRootGenerator>("root_generator")->ClearGenerator(fGenerator);
            fToolbox.Remove(fGenerator->GetName());
            fGenerator = NULL;
            LPROG(lmclog,"Removing fGenerator from KToolbox ... ");

            auto tGen = fToolbox.GetAll<Kassiopeia::KSGenerator>();
            for (unsigned i=0; i<tGen.size(); i++)
            {
                if ( (tGen[i]->IsActivated()) && (tGen[i]->GetName()!="root_generator") )
                {
                    LPROG(lmclog,"Replacing KSRoot generator with " << tGen[i]->GetName());
                    fToolbox.Get<Kassiopeia::KSRootGenerator>("root_generator")->SetGenerator(tGen[i]);
                }
            }

        }

    	return true;
    }

    bool RunPause::ExecutePostRunModification(Kassiopeia::KSRun & aRun)
    {
       	//  No interrupt has happened yet in KSRoot.  Run still in progress.
        return true;
    }


} /* namespace locust */


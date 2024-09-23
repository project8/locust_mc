/*
 * LMCTrackHold.cc
 *
 *  Created on: Apr 25, 2024
 *      Author: pslocum
 */

#include "logger.hh"
#include "LMCTrackHold.hh"


namespace locust
{

    LOGGER( lmclog, "TrackHold" );

    TrackHold::TrackHold() :
        fTrackCounter( 0 ),
        fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
        {}


    TrackHold::TrackHold( const TrackHold& aOrig ) : KSComponent(),
        fTrackCounter( 0 ),
        fInterface( aOrig.fInterface )
        {}


    TrackHold::~TrackHold()
    {
    }

    TrackHold* TrackHold::Clone() const
    {
        return new TrackHold( *this );
    }

    bool TrackHold::ConfigureByInterface()
    {

        if (fInterface->fConfigureKass)
        {
            const scarab::param_node* aParam = fInterface->fConfigureKass->GetParameters();
            if (!this->Configure( *aParam ))
            {
                LERROR(lmclog,"Error configuring TrackHold class");
                return false;
            }
        }
        else
        {
            LPROG(lmclog,"TrackHold class did not need to be configured.");
            return true;
        }
        return true;
    }

    bool TrackHold::Configure( const scarab::param_node& aParam )
    {
        return true;
    }


    bool TrackHold::ExecutePreTrackModification(Kassiopeia::KSTrack &aTrack)
    {
        fInterface->aTrack->Initialize();
        fInterface->fNewTrackStarting = true;
        double tTime = aTrack.GetInitialParticle().GetTime();

        double tPitchAngle = aTrack.GetInitialParticle().GetPolarAngleToB();
        LWARN(lmclog,"LMCTrack " << fTrackCounter << " is starting at Kass time " << tTime << " with instantaneous pitch angle " <<  tPitchAngle);
        return true;
    }

    bool TrackHold::ExecutePostTrackModification(Kassiopeia::KSTrack &aTrack)
    {
        if ( aTrack.GetTotalSteps() > 0)
        {
            fInterface->anEvent->AddTrack( fInterface->aTrack );
        }

        double tTime = aTrack.GetFinalParticle().GetTime();
        LWARN(lmclog,"LMCTrack " << fTrackCounter << " is complete at Kass time " << tTime << " with total steps " << aTrack.GetTotalSteps() );
        fTrackCounter += 1;
        return true;
    }






}  // namespace locust

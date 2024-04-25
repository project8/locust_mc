/*
 * LMCTrackHold.cc
 *
 *  Created on: Apr 25, 2024
 *      Author: pslocum
 */

#include "LMCTrackHold.hh"


namespace locust
{

TrackHold::TrackHold() {}



TrackHold::TrackHold( const TrackHold& aOrig ) {}

TrackHold::~TrackHold()
{
}

TrackHold* TrackHold::Clone() const
{
    return new TrackHold( *this );
}

bool TrackHold::ExecutePreTrackModification(Kassiopeia::KSTrack &aTrack)
{
	return true;
}

bool TrackHold::ExecutePostTrackModification(Kassiopeia::KSTrack &aTrack)
{
	return true;
}






}  // namespace locust

/*
 * LMCTrackHoldBuilder.cc
 *
 *  Created on: Apr 25, 2024
 *      Author: pslocum
 */


#ifndef LOCUST_LMCTRACKHOLD_HH_
#define LOCUST_LMCTRACKHOLD_HH_

#include "KSTrackModifier.h"
#include "KSComponentTemplate.h"

namespace locust
{

    class TrackHold :
        public Kassiopeia::KSComponentTemplate< TrackHold, Kassiopeia::KSTrackModifier >
    {
        public:
        TrackHold();
        TrackHold( const TrackHold& aOrig );
        virtual ~TrackHold();
        TrackHold* Clone() const;



        public:

        virtual bool ExecutePreTrackModification(Kassiopeia::KSTrack &aTrack);
        virtual bool ExecutePostTrackModification(Kassiopeia::KSTrack &aTrack);


};

} /* namespace locust */



#endif

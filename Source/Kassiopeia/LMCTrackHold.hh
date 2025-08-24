/*
 * LMCTrackHoldBuilder.cc
 *
 *  Created on: Apr 25, 2024
 *      Author: pslocum
 */


#ifndef LOCUST_LMCTRACKHOLD_HH_
#define LOCUST_LMCTRACKHOLD_HH_

#include "KSTrackModifier.h"
#include "KSTrack.h"

#include "KSComponentTemplate.h"
#include "LMCKassLocustInterface.hh"


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

        bool ConfigureByInterface();
        bool Configure( const scarab::param_node& aParam );




        public:

        virtual bool ExecutePreTrackModification(Kassiopeia::KSTrack &aTrack);
        virtual bool ExecutePostTrackModification(Kassiopeia::KSTrack &aTrack);

        protected:
            kl_interface_ptr_t fInterface;

        private:

        int fTrackCounter;
        bool fMottScattering;


};

} /* namespace locust */



#endif

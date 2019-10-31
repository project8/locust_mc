/*
 * LMCEvent.hh
 *
 * This class has access to both locust and ROOT libraries.  The syntax is
 * consistent with KTROOTData.hh and the instructions in
 * https://root.cern.ch/root/Using.html .  It is also mentioned in LMCEventLinkDef.hh .
 *  Created on: Dec 5, 2018
 *      Author: pslocum
 */



#ifndef LMCEVENT_HH_
#define LMCEVENT_HH_

#include "TObject.h"
#include "LMCTrack.hh"
#include <vector>

namespace locust
{

    class Event : public TObject 
    {

        public:
            Event();
            virtual ~Event();

            void AddTrack(const Track aTrack);

            int fEventID;
            double fLOFrequency;
            int fRandomSeed;

            std::vector<double> fStartFrequencies;
            std::vector<double> fTrackPowers;
            std::vector<double> fStartTimes;
            std::vector<double> fEndTimes;
            std::vector<double> fTrackLengths;
            std::vector<double> fSlopes;
            std::vector<double> fPitchAngles;

            unsigned fNTracks;

            ClassDef(Event,1)  // Root syntax.

    };

}
#endif /* LMCEVENT_HH_ */

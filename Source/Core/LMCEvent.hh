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

    class Event;
    class Event : public TObject 
    {

        public:
            Event();
            virtual ~Event();
            int EventID = -99;
            int ntracks = -99;
            double LOFrequency = -99.;
            std::vector<double> StartFrequencies;
            std::vector<double> StartTimes;
            std::vector<double> EndTimes;
            std::vector<double> TrackLengths;
            std::vector<double> Slopes;


   ClassDef(Event,1)  // Root syntax.

};

}
#endif /* LMCEVENT_HH_ */

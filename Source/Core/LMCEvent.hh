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

            int EventID;
            double LOFrequency;
            int RandomSeed;

            void AddTrack(const Track aTrack)

        private:
            std::vector<double> StartFrequencies;
            std::vector<double> TrackPower;
            std::vector<double> StartTimes;
            std::vector<double> EndTimes;
            std::vector<double> TrackLengths;
            std::vector<double> Slopes;
            std::vector<double> PitchAngles;

            unsigned ntracks;


            ClassDef(Event,1)  // Root syntax.

    };

}
#endif /* LMCEVENT_HH_ */

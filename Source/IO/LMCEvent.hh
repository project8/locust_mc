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
#include "time.h"
#include <sys/time.h>




namespace locust
{

    class Event : public TObject 
    {

        public:
            Event();
            virtual ~Event();

            bool Initialize();
            void AddTrack(const Track* aTrack);
            void AddTrack(const Track aTrack);

            long int fEventID;
            double fLOFrequency;
            long int fRandomSeed;

            std::vector<int> fTrackIDs;
            std::vector<double> fStartingEnergies_eV;
            std::vector<double> fOutputStartFrequencies;
            std::vector<double> fStartFrequencies;
            std::vector<double> fEndFrequencies;
            std::vector<double> fAvgFrequencies;
            std::vector<double> fOutputAvgFrequencies;
            std::vector<double> fAvgAxialFrequencies;
            std::vector<double> fTrackPowers;
            std::vector<double> fStartTimes;
            std::vector<double> fEndTimes;
            std::vector<double> fTrackLengths;
            std::vector<double> fSlopes;
            std::vector<double> fPitchAngles;
            std::vector<double> fRadii;
            std::vector<double> fRadialPhases;

            unsigned fNTracks;

            ClassDef(Event,1)  // Root syntax.

    };

}
#endif /* LMCEVENT_HH_ */

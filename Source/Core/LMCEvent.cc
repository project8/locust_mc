/*
 * LMCEvent.cc
 *
 * This class has access to both locust and ROOT libraries.  The syntax is
 * consistent with KTROOTData.cc and the instructions in
 * https://root.cern.ch/root/Using.html .   It is also mentioned in LMCEventLinkDef.hh .
 *
 *  Created on: Dec 5, 2018
 *      Author: pslocum
 */


#include "LMCEvent.hh"
#include <iostream>

ClassImp(locust::Event);

namespace locust
{
    Event::Event() {}
    Event::~Event() {}

    void Event::AddTrack(const Track aTrack)
    {
        fStartFrequencies.push_back( aTrack.StartFrequency );
        fTrackPowers.push_back( aTrack.TrackPower );
        fStartTimes.push_back( aTrack.StartTime );
        fTrackLengths.push_back( aTrack.TrackLength );
        fEndTimes.push_back( aTrack.EndTime );
        fSlopes.push_back( aTrack.Slope );
        fPitchAngles.push_back( aTrack.PitchAngle );

        //update size
        fNTracks = fStartFrequencies.size();
    }
}

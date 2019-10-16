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
        StartFrequencies.push_back( aTrack.StartFrequency );
        TrackPower.push_back( aTrack.TrackPower );
        StartTimes.push_back( aTrack.StartTime );
        TrackLengths.push_back( aTrack.TrackLength );
        EndTimes.push_back( aTrack.EndTime );
        Slopes.push_back( aTrack.Slope );
        PitchAngles.push_back( aTrack.PitchAngle );

        //update size
        nTracks = StartFrequencies.size();
    }

}



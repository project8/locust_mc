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
    Event::Event() :
        fNTracks( 0 )
    {
    }
    Event::~Event() {}

    void Event::AddTrack(const Track aTrack)
    {
        fStartingEnergies_eV.push_back( aTrack.StartingEnergy_eV );
    	fOutputStartFrequencies.push_back( aTrack.OutputStartFrequency );
        fStartFrequencies.push_back( aTrack.StartFrequency );
        fEndFrequencies.push_back( aTrack.EndFrequency );
        fAvgFrequencies.push_back( aTrack.AvgFrequency );
        fAvgAxialFrequencies.push_back( aTrack.AvgAxialFrequency );
        fTrackPowers.push_back( aTrack.TrackPower );
        fStartTimes.push_back( aTrack.StartTime );
        fTrackLengths.push_back( aTrack.TrackLength );
        fEndTimes.push_back( aTrack.EndTime );
        fSlopes.push_back( aTrack.Slope );
        fPitchAngles.push_back( aTrack.PitchAngle );
        fRadii.push_back( aTrack.Radius );
        fRadialPhases.push_back( aTrack.RadialPhase );

        // Update size.  And, record fLOFrequency for compatibility with previous work.  The LO frequency is
        // now also recorded in the RunParameters Tree.
        fNTracks = fStartFrequencies.size();
        fLOFrequency = aTrack.LOFrequency;
        fRandomSeed = aTrack.RandomSeed;
    }
}

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
        fNTracks( 0 ),
        fRandomSeed( 0 )
    {
    }
    Event::~Event() {}

    bool Event::Initialize()
    {
        time_t rawtime;
        struct tm * timeInfo;
        time (&rawtime);
        timeInfo = localtime (&rawtime);

        struct timeval tv;
        gettimeofday(&tv, NULL);

        int tDay = timeInfo->tm_mday;
        int tMonth = timeInfo->tm_mon + 1;
        int tYear = timeInfo->tm_year - 100;
        int tMicrosec = tv.tv_usec;

        fEventID = 1e12 + tDay*1e10 + tMonth*1e8 + tYear*1e6 + tMicrosec;
        fLOFrequency = -99.;


        fTrackIDs.resize(0);
        fStartingEnergies_eV.resize(0);
        fOutputStartFrequencies.resize(0);
        fStartFrequencies.resize(0);
        fEndFrequencies.resize(0);
        fAvgFrequencies.resize(0);
        fOutputAvgFrequencies.resize(0);
        fAvgAxialFrequencies.resize(0);
        fTrackPowers.resize(0);
        fStartTimes.resize(0);
        fTrackLengths.resize(0);
        fEndTimes.resize(0);
        fSlopes.resize(0);
        fPitchAngles.resize(0);
        fRadii.resize(0);
        fRadialPhases.resize(0);

        // Update size.
        fNTracks = 0;


        return true;
    }


    void Event::AddTrack(const Track* aTrack) // Phase III structure
    {
        fTrackIDs.push_back( aTrack->TrackID );
        fStartingEnergies_eV.push_back( aTrack->StartingEnergy_eV );
        fOutputStartFrequencies.push_back( aTrack->OutputStartFrequency );
        fStartFrequencies.push_back( aTrack->StartFrequency );
        fEndFrequencies.push_back( aTrack->EndFrequency );
        fAvgFrequencies.push_back( aTrack->AvgFrequency );
        fOutputAvgFrequencies.push_back( aTrack->OutputAvgFrequency );
        fAvgAxialFrequencies.push_back( aTrack->AvgAxialFrequency );
        fTrackPowers.push_back( aTrack->TrackPower );
        fStartTimes.push_back( aTrack->StartTime );
        fTrackLengths.push_back( aTrack->TrackLength );
        fEndTimes.push_back( aTrack->EndTime );
        fSlopes.push_back( aTrack->Slope );
        fPitchAngles.push_back( aTrack->PitchAngle );
        fRadii.push_back( aTrack->Radius );
        fRadialPhases.push_back( aTrack->RadialPhase );

        // Update size.
        fNTracks = fStartFrequencies.size();
    }

    void Event::AddTrack(const Track aTrack)  // Phase II structure
    {
        fStartingEnergies_eV.push_back( aTrack.StartingEnergy_eV );
    	fOutputStartFrequencies.push_back( aTrack.OutputStartFrequency );
        fStartFrequencies.push_back( aTrack.StartFrequency );
        fEndFrequencies.push_back( aTrack.EndFrequency );
        fAvgFrequencies.push_back( aTrack.AvgFrequency );
        fOutputAvgFrequencies.push_back( aTrack.OutputAvgFrequency );
        fAvgAxialFrequencies.push_back( aTrack.AvgAxialFrequency );
        fTrackPowers.push_back( aTrack.TrackPower );
        fStartTimes.push_back( aTrack.StartTime );
        fTrackLengths.push_back( aTrack.TrackLength );
        fEndTimes.push_back( aTrack.EndTime );
        fSlopes.push_back( aTrack.Slope );
        fPitchAngles.push_back( aTrack.PitchAngle );
        fRadii.push_back( aTrack.Radius );
        fRadialPhases.push_back( aTrack.RadialPhase );
        fNTracks = fStartFrequencies.size();
    }

}

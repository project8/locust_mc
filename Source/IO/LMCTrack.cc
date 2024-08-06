/*
 * LMCTrack.cc
 *
 *  Created on: Dec 5, 2018
 *      Author: pslocum
 */

#include "LMCTrack.hh"
ClassImp(locust::Track);


namespace locust
{

    Track::Track()
    {
    }

    Track::~Track()
    {
    }

    bool Track::Initialize()
    {
        EventID = -99;
        RandomSeed = -99.;
        StartTime = -99.;
        EndTime = -99.;
        TrackLength = -99.;
        StartingEnergy_eV = -99.;
        OutputStartFrequency = -99.;
        StartFrequency = -99.;
        EndFrequency = -99.;
        AvgFrequency = -99.;
        LOFrequency = -99.;
        TrackPower = -99.;
        Slope = -99.;
        PitchAngle = -99.;
        Radius = -99.;
        RadialPhase = -99.;

        return true;

    }


} /* namespace locust */


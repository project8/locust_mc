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

    Track::Track() :
        TrackID( -2 )
    {
    }

    Track::~Track()
    {
    }

    bool Track::Initialize()
    {
        TrackID = 0;
        StartTime = -99.;
        EndTime = -99.;
        TrackLength = -99.;
        StartingEnergy_eV = -99.;
        OutputStartFrequency = -99.;
        StartFrequency = -99.;
        EndFrequency = -99.;
        AvgFrequency = 0.;
        OutputAvgFrequency = 0.;
        TrackOutputStartFrequency = -99.;
        AvgAxialFrequency = 0.;
        TrackPower = -99.;
        Slope = -99.;
        PitchAngle = -99.;
        Radius = -99.;
        RadialPhase = -99.;
        StartGuidingCenterX = -99.;
        StartGuidingCenterY = -99.;
        StartGuidingCenterZ = -99.;

        return true;

    }


} /* namespace locust */


/*
 * LMCTrack.hh
 *
 *  Created on: Dec 5, 2018
 *      Author: pslocum
 */

#ifndef LMCTRACK_HH_
#define LMCTRACK_HH_

#include "TObject.h"


namespace locust
{
 /*!
 @class Track
 @author P. Slocum
 @brief Class to describe tracks for writing to Root file.
 @details
 Available configuration options:
 No input parameters
 */

    class Track : public TObject
    {

        public:
            Track();
            virtual ~Track();
            bool Initialize();
            int EventID = -99;
            int RandomSeed = -99.;
            double StartTime = -99.;
            double EndTime = -99.;
            double TrackLength = -99.;
            double StartingEnergy_eV = -99.;
            double OutputStartFrequency = -99.;
            double StartFrequency = -99.;
            double EndFrequency = -99.;
            double AvgFrequency = -99.;
            double AvgAxialFrequency = -99.;
            double LOFrequency = -99.;
            double TrackPower = -99.;
            double Slope = -99.;
            double PitchAngle = -99.;
            double Radius = -99.;
            double RadialPhase = -99.;

            ClassDef(Track,1)  // Root syntax.

    };

} /* namespace locust */

#endif /* LMCTRACK_HH_ */

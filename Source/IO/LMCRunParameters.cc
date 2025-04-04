/*
 * LMCRunParameters.cc
 *
 * This class has access to both locust and ROOT libraries.  The syntax is
 * consistent with KTROOTData.cc and the instructions in
 * https://root.cern.ch/root/Using.html .   It is also mentioned in LMCEventLinkDef.hh .
 *
 *  Created on: Jul 7, 2020
 *      Author: pslocum
 */


#include "LMCRunParameters.hh"
#include <iostream>

ClassImp(locust::RunParameters);

namespace locust
{
    RunParameters::RunParameters() :
        fNoise( 0. ),
        fLOfrequency( 0. ),
		fSamplingRateMHz( 0. ),
		fDecimationFactor( 0. ),
		fDataType( "simulation" ),
		fSimulationType( "rectangular-waveguide" ),
		fSimulationSubType( "none" ),
		fRunID( 0 ),
		fKassiopeiaSeed( 0 ),
		fTrackLengthSeed( 0 ),
		fTrackDelaySeed( 0 )
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

        fRunID = 1e12 + tDay*1e10 + tMonth*1e8 + tYear*1e6 + tMicrosec;

    }
    RunParameters::~RunParameters() {}

}

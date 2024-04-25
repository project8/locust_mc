/*
 * LMCTrackHoldBuilder.cc
 *
 *  Created on: Apr 25, 2024
 *      Author: pslocum
 */

#include "LMCTrackHoldBuilder.hh"

#include "KSRootBuilder.h"

template< >
locust::TrackHoldBuilder::~KComplexElement()
{
}

namespace locust
{

    STATICINT SLMCTrackHoldStructure =
            locust::TrackHoldBuilder::Attribute< std::string >( "name" );

    STATICINT sLMCTrackHold =
            katrin::KSRootBuilder::ComplexElement< locust::TrackHold >( "mod_track_hold" );

} /* namespace locust */


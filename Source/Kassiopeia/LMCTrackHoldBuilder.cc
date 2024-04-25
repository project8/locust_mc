/*
 * LMCEventHoldBuilder.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCEventHoldBuilder.hh"

#include "KSRootBuilder.h"

template< >
locust::EventHoldBuilder::~KComplexElement()
{
}

namespace locust
{

    STATICINT SLMCEventHoldStructure =
            locust::EventHoldBuilder::Attribute< std::string >( "name" );

    STATICINT sLMCEventHold =
            katrin::KSRootBuilder::ComplexElement< locust::EventHold >( "mod_event_hold" );

} /* namespace locust */


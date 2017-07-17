/*
 * LMCEventHoldBuilder.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCEventHoldBuilder.hh"

#include "KSRootBuilder.h"

namespace katrin
{

    template< >
    EventHoldBuilder::~KComplexElement()
    {
    }

    STATICINT SLMCEventHoldStructure =
            EventHoldBuilder::Attribute< std::string >( "name" )+
            EventHoldBuilder::Attribute< bool >( "wait_before_event" )+
            EventHoldBuilder::Attribute< bool >( "wait_after_event" );

    STATICINT sLMCEventHold =
            KSRootBuilder::ComplexElement< locust::EventHold >( "mod_event_hold" );

} /* namespace katrin */

/*
 * LMCEventHoldBuilder.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCEventHoldBuilder.hh"

using namespace katrin;

namespace locust
{

    template< >
    EventHoldBuilder::~KComplexElement()
    {
    }

    STATICINT SLMCEventHoldStructure =
            EventHoldBuilder::Attribute< string >( "name" )+
            EventHoldBuilder::Attribute< bool >( "wait_before_event" )+
            EventHoldBuilder::Attribute< string >( "wait_after_event" );

    STATICINT sLMCEventHold =
            KSRootBuilder::ComplexElement< EventHold >( "mod_event_hold" );

} /* namespace locust */

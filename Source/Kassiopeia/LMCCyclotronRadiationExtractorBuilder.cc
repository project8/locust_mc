/*
 * LMCCyclotronRadiationExtractorBuilder.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCCyclotronRadiationExtractorBuilder.hh"

using namespace katrin;

namespace locust
{

    template< >
    CyclotronRadiationExtractorBuilder::~KComplexElement()
    {
    }

    STATICINT SLMCEventHoldStructure =
            CyclotronRadiationExtractorBuilder::Attribute< string >( "name" )/*+
            CyclotronRadiationExtractorBuilder::Attribute< bool >( "wait_before_event" )+
            CyclotronRadiationExtractorBuilder::Attribute< string >( "wait_after_event" )*/;

    STATICINT sLMCEventHold =
            KSRootBuilder::ComplexElement< CyclotronRadiationExtractor >( "cycl_rad_extr" );

} /* namespace locust */

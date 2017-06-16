/*
 * LMCCyclotronRadiationExtractorBuilder.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCCyclotronRadiationExtractorBuilder.hh"

#include "KSRootBuilder.h"

namespace katrin
{

    template< >
    CyclotronRadiationExtractorBuilder::~KComplexElement()
    {
    }

    STATICINT SLMCCyclRadExtrStructure =
            CyclotronRadiationExtractorBuilder::Attribute< string >( "name" )/*+
            CyclotronRadiationExtractorBuilder::Attribute< bool >( "wait_before_event" )+
            CyclotronRadiationExtractorBuilder::Attribute< string >( "wait_after_event" )*/;

    STATICINT sLMCCyclotronRadiationExtractor =
            KSRootBuilder::ComplexElement< locust::CyclotronRadiationExtractor >( "cycl_rad_extr" );

} /* namespace katrin */

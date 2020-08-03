/*
 * LMCRunPauseBuilder.cc
 *
 *  Created on: Aug 28, 2019
 *      Author: N.S. Oblath
 */

#include "LMCRunPauseBuilder.hh"

#include "KSRootBuilder.h"

template< >
locust::RunPauseBuilder::~KComplexElement()
{
}

namespace locust
{

    STATICINT SLMCRunPauseStructure =
            locust::RunPauseBuilder::Attribute< std::string >( "name" );

    STATICINT sLMCRunPause =
            katrin::KSRootBuilder::ComplexElement< locust::RunPause >( "mod_run_pause" );

} /* namespace locust */


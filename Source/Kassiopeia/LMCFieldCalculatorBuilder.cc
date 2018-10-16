/*
 * LMCFieldCalculatorBuilder.cc
 *
 *  Created on: Oct 4, 2018
 *      Author: pslocum
 */

#include "LMCFieldCalculatorBuilder.hh"

#include "KSRootBuilder.h"
using namespace std;

using namespace Kassiopeia;
namespace katrin
{

    template< >
    FieldCalculatorBuilder::~KComplexElement()
    {
    }

    STATICINT SLMCFieldCalcStructure =
      FieldCalculatorBuilder::Attribute< std::string >( "name" );

    STATICINT sLMCFieldCalculator =
            KSRootBuilder::ComplexElement< locust::FieldCalculator >( "field_calc" );

} /* namespace katrin */


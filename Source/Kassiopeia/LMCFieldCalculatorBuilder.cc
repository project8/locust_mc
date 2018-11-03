#include "LMCFieldCalculatorBuilder.hh"
#include "KSRootBuilder.h"

using namespace Kassiopeia;
using namespace std;

namespace katrin
{

    template< >
    FieldCalculatorBuilder::~KComplexElement()
    {
    }

    STATICINT sFieldCalculator =
      KSRootBuilder::ComplexElement< locust::FieldCalculator >( "field_calculator" );

    STATICINT sFieldCalculatorStructure =
        FieldCalculatorBuilder::Attribute< string >( "name" ) +
        FieldCalculatorBuilder::Attribute< string >( "add_space_interaction" );

}

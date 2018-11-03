#include "LMCCyclotronRadiationExtractorBuilder.hh"
#include "KSRootBuilder.h"

using namespace Kassiopeia;
using namespace std;

namespace katrin
{
    template< >
    CyclotronRadiationExtractorBuilder::~KComplexElement()
    {
    }

    STATICINT sLMCCyclotronRadiationExtractor =
      KSRootBuilder::ComplexElement< CyclotronRadiationExtractor >( "cycl_rad_extr" );

    STATICINT SLMCCyclRadExtrStructure =
            CyclotronRadiationExtractorBuilder::Attribute< string >( "name" ) +
            CyclotronRadiationExtractorBuilder::Attribute< int >( "P8Phase" ) +
            CyclotronRadiationExtractorBuilder::Attribute< string >( "add_modifier" );
}

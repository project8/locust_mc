#include "LMCRootBuilder.hh"
#include "KElementProcessor.hh"

using namespace Kassiopeia;
namespace katrin
{

    template< >
    LMCRootBuilder::~KComplexElement()
    {
    }

    static int sLMCRoot =
        KElementProcessor::ComplexElement< KSRoot >( "lmc_kassiopeia" );

}

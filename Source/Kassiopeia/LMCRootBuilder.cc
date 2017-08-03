#include "LMCRootBuilder.hh"
#include "KElementProcessor.hh"
#include "KRoot.h"


using namespace Kassiopeia;
using namespace std;

namespace katrin
{

    template< >
    LMCRootBuilder::~KComplexElement()
    {
    }

    STATICINT sLMCRoot =
//        KElementProcessor::ComplexElement< KSRoot >( "lmc_kassiopeia" );
    KRootBuilder::ComplexElement< KSRoot >( "lmc_kassiopeia" );



    STATICINT sLMCRootCompat =
        KElementProcessor::ComplexElement< KSRoot >( "lmc_kassiopeia" );


}

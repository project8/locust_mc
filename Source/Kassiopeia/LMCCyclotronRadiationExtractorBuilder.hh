
#ifndef LOCUST_LMCCYCLOTRONRADIATIONEXTRACTORBUILDER_HH_
#define LOCUST_LMCCYCLOTRONRADIATIONEXTRACTORBUILDER_HH_



#include "KComplexElement.hh"
#include "LMCCyclotronRadiationExtractor.hh"
#include "KToolbox.h"

using namespace Kassiopeia;
using namespace locust;
namespace katrin
{
typedef KComplexElement< CyclotronRadiationExtractor > CyclotronRadiationExtractorBuilder;

    template< >
    inline bool CyclotronRadiationExtractorBuilder::AddAttribute(KContainer *aContainer)
    {
        if( aContainer->GetName() == "name" )
        {
            aContainer->CopyTo( fObject, &KNamed::SetName );
            return true;
        }
        if( aContainer->GetName() == "P8Phase" )
        {
	    aContainer->CopyTo( fObject, &CyclotronRadiationExtractor::SetP8Phase );
            return true;
        }
        if( aContainer->GetName() == "add_modifier" )
        {
            fObject->AddModifier( KToolbox::GetInstance().Get< KSStepModifier >( aContainer->AsReference< std::string >() ) );
            return true;
        }
        return false;
    }
}

#endif //Kassiopeia_KSRootStepModifierBuilder_h_

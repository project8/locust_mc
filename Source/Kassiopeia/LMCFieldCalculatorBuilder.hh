#ifndef LOCUST_LMCFIELDCALCULATORBUILDER_HH_
#define LOCUST_LMCFIELDCALCULATORBUILDER_HH_

#include "KComplexElement.hh"
#include "LMCFieldCalculator.hh"
#include "KToolbox.h"

using namespace Kassiopeia;
namespace katrin
{

  typedef KComplexElement< locust::FieldCalculator > FieldCalculatorBuilder;

    template< >
    inline bool FieldCalculatorBuilder::AddAttribute( KContainer* aContainer )
    {
        if( aContainer->GetName() == "name" )
        {
            aContainer->CopyTo( fObject, &KNamed::SetName );
            return true;
        }
        if( aContainer->GetName() == "add_space_interaction" )
        {
            fObject->AddSpaceInteraction( KToolbox::GetInstance().Get< KSSpaceInteraction >( aContainer->AsReference< std::string >() ) );
            return true;
        }
        return false;
    }

}
#endif

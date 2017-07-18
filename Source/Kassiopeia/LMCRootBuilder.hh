#ifndef LMC_KSRootBuilder_h_
#define LMC_KSRootBuilder_h_

#include "KComplexElement.hh"
#include "KSRoot.h"
#include "KSSimulation.h"
#include "KToolbox.h"

using namespace Kassiopeia;
namespace katrin
{

    typedef KComplexElement< KSRoot > LMCRootBuilder;

    template< >
    inline bool LMCRootBuilder::Begin()
    {
        fObject = new KSRoot();
        return true;
    }

    template< >
    inline bool LMCRootBuilder::AddElement( KContainer* aContainer )
    {
        if( aContainer->Is< KSSimulation >() )
        {
            std::cout << "### running object called <" << fObject->GetName() << ">" << std::endl;
            aContainer->ReleaseTo( fObject, &KSRoot::Execute );
	    //            printf("running the object check 1\n");
            return true;
        }
        if( aContainer->Is< KSObject >() )
        {
            std::cout << "### adding object called <" << fObject->GetName() << ">" << std::endl;
            KToolbox::GetInstance().AddContainer(*aContainer);
            return true;
        }
        return false;
    }

    template< >
    inline bool LMCRootBuilder::End()
    {
        delete fObject;
        return true;
    }

}

#endif

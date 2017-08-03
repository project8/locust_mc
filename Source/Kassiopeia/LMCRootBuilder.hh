#ifndef LMC_KSRootBuilder_h_
#define LMC_KSRootBuilder_h_

#include "KComplexElement.hh"
#include "KSRoot.h"
#include "KSSimulation.h"
#include "KToolbox.h"

#include "LMCCyclotronRadiationExtractor.hh"

#include "KEMToolbox.hh"
#include "KElectricField.hh"
#include "KSElectricField.h"
#include "KSElectricKEMField.h"
#include "KMagneticField.hh"
#include "KSMagneticField.h"
#include "KSMagneticKEMField.h"
#include "KElectrostaticPotentialmapBuilder.hh"


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
            aContainer->ReleaseTo( fObject, &KSRoot::Execute );
            return true;
        }
        if( aContainer->Is< KSObject >() )
        {
<<<<<<< HEAD
            std::cout << "### adding object called <" << fObject->GetName() << ">" << std::endl;
            KToolbox::GetInstance().AddContainer(*aContainer);
=======
            KToolbox::GetInstance().AddContainer(*aContainer);
            return true;

        }

        // legacy support for old field bindings in the <kassiopeia> tag
        if( aContainer->Is< KEMField::KElectricField >() ||
            aContainer->Is< KEMField::KElectrostaticPotentialmapCalculator >() )
        {
            KSElectricKEMField* tField = new KSElectricKEMField();
            tField->SetName(aContainer->GetName());
            aContainer->ReleaseTo(tField, &KSElectricKEMField::SetElectricField );
            KToolbox::GetInstance().Add(tField,tField->GetName());
            return true;
        }
        if( aContainer->Is< KEMField::KMagneticField >() )
        {
            KSMagneticKEMField* tField = new KSMagneticKEMField();
            tField->SetName(aContainer->GetName());
            aContainer->ReleaseTo(tField, &KSMagneticKEMField::SetMagneticField );
            KToolbox::GetInstance().Add(tField,tField->GetName());
>>>>>>> my-temporary-work
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

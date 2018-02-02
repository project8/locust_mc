/*
 * LMCCyclotronRadiationExtractorBuilder.hh
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#ifndef LOCUST_LMCCYCLOTRONRADIATIONEXTRACTORBUILDER_HH_
#define LOCUST_LMCCYCLOTRONRADIATIONEXTRACTORBUILDER_HH_

#include "KComplexElement.hh"
#include "LMCCyclotronRadiationExtractor.hh"
#include "KToolbox.h"

using namespace locust;
namespace katrin
{
    /*!
     @class 
     @author N. S. Oblath

     @brief bindings for CyclotronExtractor for the locust xml file

     @details

     Available configuration options:
     - "name": string -- name of xml object
     - "P8Phase": int -- phase of P8 experiment. Currently 1, 2, 3 are supported
     - "set_trajectrory": string -- used for interpolation
    */

    typedef KComplexElement< locust::CyclotronRadiationExtractor > CyclotronRadiationExtractorBuilder;

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
        if( aContainer->GetName() == "set_trajectory" )
        {
            fObject->SetTrajectory( KToolbox::GetInstance().Get< Kassiopeia::KSTrajectory >( aContainer->AsReference< std::string >() ) );
            return true;
        }

        /*
        if( aContainer->GetName() == "wait_before_event" )
        {
            fObject->SetWaitBeforeEvent( aContainer->AsReference< bool >() );
            return true;
        }
        if( aContainer->GetName() == "wait_after_event" )
        {
            fObject->SetWaitAfterEvent( aContainer->AsReference< bool >() );
            return true;
        }
        */
        return false;
    }

} /* namespace katrin */

#endif /* LOCUST_LMCCYCLOTRONRADIATIONEXTRACTORBUILDER_HH_ */


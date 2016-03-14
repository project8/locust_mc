/*
 * LMCEventHoldBuilder.hh
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#ifndef LOCUST_LMCEVENTHOLDBUILDER_HH_
#define LOCUST_LMCEVENTHOLDBUILDER_HH_

#include "KComplexElement.hh"

#include "LMCEventHold.hh"

namespace locust
{

    typedef katrin::KComplexElement< EventHold > EventHoldBuilder;

    template< >
    inline bool EventHoldBuilder::AddAttribute(KContainer *aContainer)
    {
        if( aContainer->GetName() == "name" )
        {
            aContainer->CopyTo( fObject, &KNamed::SetName );
            return true;
        }
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
        return false;
    }

} /* namespace locust */

#endif /* LOCUST_LMCEVENTHOLDBUILDER_HH_ */

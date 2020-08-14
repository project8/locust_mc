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

    typedef katrin::KComplexElement< locust::EventHold > EventHoldBuilder;

} /* namespace locust */

template< >
inline bool locust::EventHoldBuilder::AddAttribute(KContainer *aContainer)
{
    if( aContainer->GetName() == "name" )
    {
        aContainer->CopyTo( fObject, &KNamed::SetName );
        return true;
    }
    return false;
}


#endif /* LOCUST_LMCEVENTHOLDBUILDER_HH_ */

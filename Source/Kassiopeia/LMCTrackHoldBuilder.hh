/*
 * LMCTrackHoldBuilder.hh
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#ifndef LOCUST_LMCTRACKHOLDBUILDER_HH_
#define LOCUST_LMCTRACKHOLDBUILDER_HH_

#include "KComplexElement.hh"

#include "LMCTrackHold.hh"

namespace locust
{

    typedef katrin::KComplexElement< locust::TrackHold > TrackHoldBuilder;

} /* namespace locust */

template< >
inline bool locust::TrackHoldBuilder::AddAttribute(KContainer *aContainer)
{
    if( aContainer->GetName() == "name" )
    {
        aContainer->CopyTo( fObject, &KNamed::SetName );
        return true;
    }
    return false;
}


#endif /* LOCUST_LMCTRACKHOLDBUILDER_HH_ */

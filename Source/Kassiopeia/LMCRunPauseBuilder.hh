/*
 * LMCRunPauseBuilder.hh
 *
 *  Created on: Aug 28, 2019
 *      Author: N.S. Oblath
 */

#ifndef LOCUST_LMCRUNPAUSEBUILDER_HH_
#define LOCUST_LMCRUNPAUSEBUILDER_HH_

#include "KComplexElement.hh"

#include "LMCRunPause.hh"

namespace locust
{

    typedef katrin::KComplexElement< locust::RunPause > RunPauseBuilder;

} /* namespace locust */

template< >
inline bool locust::RunPauseBuilder::AddAttribute(KContainer *aContainer)
{
    if( aContainer->GetName() == "name" )
    {
        aContainer->CopyTo( fObject, &KNamed::SetName );
        return true;
    }
    return false;
}

#endif /* LOCUST_LMCRUNPAUSEBUILDER_HH_ */

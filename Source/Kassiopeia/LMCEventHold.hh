/*
 * LMCEventHold.hh
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#ifndef LOCUST_LMCEVENTHOLD_HH_
#define LOCUST_LMCEVENTHOLD_HH_

#include "KSEventModifier.h"
#include "KSComponentTemplate.h"

#include "LMCKassLocustInterface.hh"

namespace locust
{

    class EventHold :
            public Kassiopeia::KSComponentTemplate< EventHold, Kassiopeia::KSEventModifier >
    {
        public:
            EventHold();
            EventHold( const EventHold& aOrig );
            virtual ~EventHold();

            EventHold* Clone() const;

        public:

            virtual bool ExecutePreEventModification(Kassiopeia::KSEvent &anEvent);
            virtual bool ExecutePostEventModification(Kassiopeia::KSEvent &anEvent);

        protected:
            kl_interface_ptr_t fInterface;

    };

} /* namespace locust */

#endif /* LOCUST_LMCEVENTHOLD_HH_ */


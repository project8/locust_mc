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

#ifdef ROOT_FOUND
    #include "LMCRootTreeWriter.hh"
    #include "LMCEvent.hh"
#endif


namespace locust
{

    class EventHold :
            public Kassiopeia::KSComponentTemplate< EventHold, Kassiopeia::KSEventModifier >
    {
        public:
            EventHold();
            EventHold( const EventHold& aOrig );
            bool OpenEvent();
            bool WriteEvent();
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


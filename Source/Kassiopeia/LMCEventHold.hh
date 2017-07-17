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

#include "KField.h"

#include <condition_variable>
#include <mutex>

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

            K_SET_GET( bool, WaitBeforeEvent );
            K_SET_GET( bool, WaitAfterEvent );

        public:

            virtual bool ExecutePreEventModification( Kassiopeia::KSEvent& );
            virtual bool ExecutePostEventModification( Kassiopeia::KSEvent& );

        public:
            void WakeBeforeEvent();
            void WakeAfterEvent();

        private:
            std::mutex fMutex;
            std::condition_variable fPreEventCondition;
            std::condition_variable fPostEventCondition;

        private:
            void InitializeComponent();
            void DeinitializeComponent();

        protected:
            virtual void PullDeupdateComponent();
            virtual void PushDeupdateComponent();


    };

} /* namespace locust */

#endif /* LOCUST_LMCEVENTHOLD_HH_ */

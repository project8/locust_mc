/*
 * LMCEventHold.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCEventHold.hh"

namespace locust
{

    EventHold::EventHold() :
            fWaitBeforeEvent( false ),
            fWaitAfterEvent( false ),
            fMutex(),
            fPreEventCondition(),
            fPostEventCondition()
    {
    }

    EventHold::EventHold( const EventHold& aOrig ) :
            fWaitBeforeEvent( aOrig.fWaitBeforeEvent ),
            fWaitAfterEvent( aOrig.fWaitAfterEvent ),
            fMutex(),
            fPreEventCondition(),
            fPostEventCondition()
    {
    }

    EventHold::~EventHold()
    {
    }

    EventHold* EventHold::Clone() const
    {
        return new EventHold( *this );
    }


    bool EventHold::ExecutePreEventModification(Kassiopeia::KSEvent &anEvent)
    {

    	printf("check 1\n"); getchar();
        if( fWaitBeforeEvent )
        {
            std::unique_lock< std::mutex >tLock( fMutex );
            fPreEventCondition.wait( tLock );
            return true;
        }
        return false;
    }

    bool EventHold::ExecutePostEventModification(Kassiopeia::KSEvent &anEvent)
    {
        if( fWaitAfterEvent )
        {
            std::unique_lock< std::mutex >tLock( fMutex );
            fPostEventCondition.wait( tLock );
            return true;
        }
        return false;
    }

    void EventHold::WakeBeforeEvent()
    {
        fPreEventCondition.notify_one();
        return;
    }

    void EventHold::WakeAfterEvent()
    {
        fPostEventCondition.notify_one();
        return;
    }

    void EventHold::InitializeComponent()
    {
    }

    void EventHold::DeinitializeComponent()
    {
    }

    void EventHold::PullDeupdateComponent()
    {
    }
    void EventHold::PushDeupdateComponent()
    {
    }


} /* namespace locust */

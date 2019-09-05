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
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    EventHold::EventHold( const EventHold& aOrig ) : KSComponent(),
            fInterface( aOrig.fInterface )
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
        printf("Kass is waiting for event trigger.\n");

        fInterface->fKassEventReady = true;
        fInterface->fFalseStartKassiopeia = false;
        fInterface->fDigitizerCondition.notify_one();  // unlock if still locked.
        if( fInterface->fWaitBeforeEvent )  // true by default
        {
            fInterface->fKassReadyCondition.notify_one();
            std::cout << "going to wait on the pre-event condition now" << std::endl;
            std::unique_lock< std::mutex >tLock( fInterface->fMutex );
            fInterface->fPreEventCondition.wait( tLock );
            fInterface->fKassEventReady = false;
            fInterface->fEventInProgress = true; // possibly redundant.
            fInterface->fTOld = 0.;  // reset time on digitizer.
        }

        printf("Kass got the event trigger\n");
        return true;

    }

    bool EventHold::ExecutePostEventModification(Kassiopeia::KSEvent &anEvent)
    {
        fInterface->fEventInProgress = false;
        fInterface->fDigitizerCondition.notify_one();  // unlock
        printf("Kass is waking after event\n");
        return true;
    }

} /* namespace locust */


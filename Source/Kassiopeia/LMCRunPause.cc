/*
 * LMCRunPause.cc
 *
 *  Created on: Jul 31, 2019
 *      Author: N.S. Oblath
 */

#include "LMCRunPause.hh"

#include "KSRun.h"

#include <csignal>


namespace locust
{

    RunPause::RunPause() :
            fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }

    RunPause::RunPause( const RunPause& aCopy ) :
            KSComponent(),
            fInterface( aCopy.fInterface )
    {
    }

    RunPause::~RunPause()
    {
    }

    RunPause* RunPause::Clone() const
    {
        return new RunPause( *this );
    }


    bool RunPause::ExecutePreRunModification(Kassiopeia::KSRun &)
    {
    	/*
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
        if (! fInterface->fRunInProgress)
        {
            printf("Raising sigint to cancel Kassiopeia");
            raise(SIGINT);
            return false;
        }
        */
        return true;
    }

    bool RunPause::ExecutePostRunModification(Kassiopeia::KSRun &)
    {
        std::cout << "executing post-run mod now" << std::endl;
        fInterface->fRunInProgress = false;
        fInterface->fKassReadyCondition.notify_one();
        return true;
    }


} /* namespace locust */


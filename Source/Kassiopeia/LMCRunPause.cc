/*
 * LMCRunPause.cc
 *
 *  Created on: Jul 31, 2019
 *      Author: N.S. Oblath
 */

#include "LMCRunPause.hh"

#include <csignal>


namespace locust
{

    RunPause::RunPause( kl_interface_ptr_t aInterface ) :
            fInterface( aInterface )
    {
    }

    RunPause::RunPause( const RunPause& ) :
            KSComponent()
    {
    }

    RunPause::~RunPause()
    {
    }

    RunPause* RunPause::Clone() const
    {
        return new RunPause( *this );
    }


    bool RunPause::ExecutePreEventModification(Kassiopeia::KSRun &aRun)
    {
        printf("Kass is waiting for event trigger.\n");


        fInterface->fKassEventReady = true;
        fInterface->fFalseStartKassiopeia = false;
        fInterface->fDigitizerCondition.notify_one();  // unlock if still locked.
        if( fInterface->fWaitBeforeEvent )  // true by default
        {
            fInterface->fKassReadyCondition.notify_one();
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
        return true;
    }

    bool RunPause::ExecutePostEventModification(Kassiopeia::KSRun &aRun)
    {
        WakeAfterEvent(fSimulation->GetEvents(), fRun->GetTotalEvents());
        return true;
    }

    void RunPause::WakeAfterEvent(unsigned TotalEvents, unsigned EventsSoFar)
    {
        fInterface->fEventInProgress = false;
        if( TotalEvents == EventsSoFar )
        {
            fInterface->fRunInProgress = false;
            fInterface->fKassReadyCondition.notify_one();
        }
        fInterface->fDigitizerCondition.notify_one();  // unlock
        printf("Kass is waking after event\n");
        return;
    }



    void RunPause::InitializeComponent()
    {
    }

    void RunPause::DeinitializeComponent()
    {
    }

    void RunPause::PullDeupdateComponent()
    {
    }
    void RunPause::PushDeupdateComponent()
    {
    }


} /* namespace locust */


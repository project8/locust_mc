/*
 * LMCRunPause.cc
 *
 *  Created on: Jul 31, 2019
 *      Author: N.S. Oblath
 */

#include "LMCRunPause.hh"

#include "LMCGlobalsDeclaration.hh"
#include "LMCGlobalsDefinition.hh"

#include <csignal>


namespace locust
{

    RunPause::RunPause()
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


        fKassEventReady = true;
        fFalseStartKassiopeia = false;
        fDigitizerCondition.notify_one();  // unlock if still locked.
        if( fWaitBeforeEvent )  // true by default
        {
            fKassReadyCondition.notify_one();
            std::unique_lock< std::mutex >tLock( fMutex );
            fPreEventCondition.wait( tLock );
            fKassEventReady = false;
            fEventInProgress = true; // possibly redundant.
            t_old = 0.;  // reset time on digitizer.
        }

        printf("Kass got the event trigger\n");
        if (! fRunInProgress)
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
        fEventInProgress = false;
        if( TotalEvents == EventsSoFar )
        {
            fRunInProgress = false;
            fKassReadyCondition.notify_one();
        }
        fDigitizerCondition.notify_one();  // unlock
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


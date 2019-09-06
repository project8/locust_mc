/*
 * LMCRunPause.cc
 *
 *  Created on: Jul 31, 2019
 *      Author: N.S. Oblath
 */

#include "LMCRunPause.hh"

#include "KSRun.h"


#include "KToolbox.h"

#include <csignal>

using namespace katrin;
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
            fInterface->fTOld = 0.;  // reset time on digitizer.
        }

        printf("Kass got the event trigger\n");
        if (! fInterface->fRunInProgress)
        {
//         printf("Raising sigint to cancel Kassiopeia");
//         raise(SIGINT);
         return false;
        }

        return true;
    }

    bool RunPause::ExecutePostRunModification(Kassiopeia::KSRun &)
    {
        Kassiopeia::KSRun* fRun = KToolbox::GetInstance().Get<Kassiopeia::KSRun>("run");

// still need to access fSimulation->GetEvents() for this to work.
//        if ( fRun->GetTotalEvents() < fSimulation->GetEvents())
        if ( fRun->GetTotalEvents() < 1)
        {
        	return true;
        }
        else
        {
        	fInterface->fRunInProgress = false;
        }


        return true;
    }


} /* namespace locust */


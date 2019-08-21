/*
 * LMCKassLocustInterface.cc
 *
 *  Created on: Aug 21, 2019
 *      Author: N.S. Oblath
 */

#include "LMCKassLocustInterface.hh"

namespace locust
{

    KassLocustInterface::KassLocustInterface() :
            fTOld( -99. ),
            fKassTimeStep( 0. ),
            fParticleHistory(),
            fWaitBeforeEvent( true ),
            fWaitAfterEvent( true ),
            fKassEventReady( false ),
            fEventInProgress( false ),
            fRunInProgress( false ),
            fPreEventInProgress( false ),
            fFalseStartKassiopeia( true ),
            fDoneWithSignalGeneration( false ),
            fMutex(),
            fKassReadyMutex(),
            fMutexDigitizer(),
            fPreEventCondition(),
            fPostEventCondition(),
            fDigitizerCondition(),
            fKassReadyCondition()
    {
    }

} /* namespace locust */

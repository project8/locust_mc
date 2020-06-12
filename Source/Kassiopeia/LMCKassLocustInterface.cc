/*
 * LMCKassLocustInterface.cc
 *
 *  Created on: Aug 21, 2019
 *      Author: N.S. Oblath
 */

#include "LMCKassLocustInterface.hh"
#include <exception>


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
            fPreEventInProgress( false ),
            fFalseStartKassiopeia( true ),
            fDoneWithSignalGeneration( false ),
            fMutex(),
            fKassReadyMutex(),
            fMutexDigitizer(),
            fPreEventCondition(),
            fPostEventCondition(),
            fDigitizerCondition(),
            fKassReadyCondition(),
            fProject8Phase( 0 ),
            fCENTER_TO_SHORT( 0.05 ),
            fCENTER_TO_ANTENNA( 0.05 )
    {}

    KLInterfaceBootstrapper::KLInterfaceBootstrapper() :
            fInterface()
    {}

    KLInterfaceBootstrapper::~KLInterfaceBootstrapper()
    {}

    kl_interface_ptr_t KLInterfaceBootstrapper::GetInterface() const
    {
        if( ! fInterface )
        {
            throw std::runtime_error( "Interface has not been set" );
        }
        return fInterface;
    }

    void KLInterfaceBootstrapper::SetInterface( kl_interface_ptr_t aInterface )
    {
        fInterface = aInterface;
        return;
    }

} /* namespace locust */

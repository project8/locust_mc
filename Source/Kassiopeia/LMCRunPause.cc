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


        return true;
    }

    bool RunPause::ExecutePostRunModification(Kassiopeia::KSRun &)
    {

    	// for some reason this is being executed after every event.

        return true;
    }


} /* namespace locust */


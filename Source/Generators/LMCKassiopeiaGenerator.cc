/*
 * LMCKassiopeiaGenerator.cc
 *
 *  Created on: June 29, 2015
 *      Author: plslocum after nsoblath
 */

#include "LMCKassiopeiaGenerator.hh"
#include "../Core/LMCLogger.hh"

#include "KSRoot.h"
#include "KSSimulation.h"
#include "KMessage.h"
#include "KSRunMessage.h"
#include "KSEventMessage.h"
#include "KSTrackMessage.h"
#include "KSStepMessage.h"


using std::string;

namespace locust
{
    LMCLOGGER( lmclog, "KassiopeiaGenerator" );

    MT_REGISTER_GENERATOR(KassiopeiaGenerator, "kassiopeia");

    KassiopeiaGenerator::KassiopeiaGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &KassiopeiaGenerator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::kTime;
    }

    KassiopeiaGenerator::~KassiopeiaGenerator()
    {
    }

    bool KassiopeiaGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;
        if( aParam->Has( "domain" ) )
        {
            string domain = aParam->GetValue( "domain" );
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LMCDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LMCERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
                return false;
            }
        }

        return true;
    }

    void KassiopeiaGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    Signal::State KassiopeiaGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void KassiopeiaGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &KassiopeiaGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &KassiopeiaGenerator::DoGenerateFreq;
        }
        else
        {
            LMCWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }

    bool KassiopeiaGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool KassiopeiaGenerator::DoGenerateTime( Signal* aSignal ) const
    {

        Kassiopeia::mainmsg( katrin::eNormal ) << "Hello world from Kassiopeia, inside Locust. " << katrin::eom;
        getchar();
//        Kassiopeia::KSSimulation* aSimulation;
//        printf("Here I am!\n");  getchar();
        return true;
    }

    bool KassiopeiaGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
         return true;
    }

} /* namespace locust */

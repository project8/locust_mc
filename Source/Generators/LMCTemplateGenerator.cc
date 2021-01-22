/*
 * LMC[name]Generator.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#include "LMC[name]Generator.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "[name]Generator" );

    MT_REGISTER_GENERATOR([name]Generator, "config-name");

    [name]Generator::[name]Generator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &[name]Generator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::k[domain];
    }

    [name]Generator::~[name]Generator()
    {
    }

    bool [name]Generator::Configure( const scarab::param_node& aParam )
    {
        return true;
    }

    void [name]Generator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    Signal::State [name]Generator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void [name]Generator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &[name]Generator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &[name]Generator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    bool [name]Generator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool [name]Generator::DoGenerateTime( Signal* aSignal )
    {
//        aSignal->SignalTimeComplex()[index][0] +=  [...]
//        aSignal->SignalTimeComplex()[index][1] +=  [...]
        return true;
    }

    bool [name]Generator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }


} /* namespace locust */

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
            Generator( aName )
    {
        fRequiredSignalState = Signal::k[domain];
    }

    [name]Generator::~[name]Generator()
    {
    }

    bool [name]Generator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        return true;
    }

    void [name]Generator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    bool [name]Generator::DoGenerate( Signal* aSignal ) const
    {
        for( unsigned index = 0; index < aSignal->[domain]Size(); ++index )
        {
            aSignal->Signal[domain]S( index ) += ???;
        }
        return true;
    }

} /* namespace locust */

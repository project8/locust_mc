/*
 * LMCGeneratorToolbox.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#include "LMCGeneratorToolbox.hh"

#include "LMCFactory.hh"
#include "LMCGenerator.hh"
#include "LMCLogger.hh"
#include "LMCParam.hh"

namespace locust
{
    LMCLOGGER( lmclog, "GeneratorToolbox" );

    GeneratorToolbox::GeneratorToolbox() :
            fFirstGenerator( NULL ),
            fRNG()
    {
        fRNG.Reseed();
    }

    GeneratorToolbox::~GeneratorToolbox()
    {
        Generator* nextGenerator = fFirstGenerator;
        while( nextGenerator != NULL )
        {
            Generator* thisGenerator = nextGenerator;
            LMCDEBUG( lmclog, "Cleaning up " << thisGenerator->GetName() );
            nextGenerator = thisGenerator->GetNextGenerator();
            delete thisGenerator;
        }
    }

    void GeneratorToolbox::Configure( const ParamNode* aNode )
    {
        if( aNode == NULL ) return;

        LMCINFO( lmclog, "Creating generators" );

        Factory< Generator >* genFactory = Factory< Generator >::GetInstance();

        const ParamArray* generatorList = aNode->array_at( "generator-list" );
        if( generatorList == NULL )
        {
            LMCERROR( lmclog, "No generator list was found" );
            return;
        }

        Generator* lastGenerator = NULL;
        for( ParamArray::const_iterator it = generatorList->begin(); it != generatorList->end(); ++it )
        {
            if( ! (*it)->is_value() )
            {
                LMCERROR( lmclog, "Non-value-type array element found in generator-list" );
                continue;
            }
            Generator* newGenerator = genFactory->Create( (*it)->as_value().get() );
            if( newGenerator == NULL )
            {
                LMCERROR( lmclog, "Unrecognized generator name: " << (*it)->as_value().get() );
                continue;
            }
            if( lastGenerator == NULL )
            {
                fFirstGenerator = newGenerator;
            }
            else
            {
                lastGenerator->SetNextGenerator( newGenerator );
            }
        }

        LMCINFO( lmclog, "Configuring generators" );

        Generator* nextGenerator = fFirstGenerator;
        while( nextGenerator != NULL )
        {
            LMCINFO( lmclog, "Configuring " << nextGenerator->GetName() );
            nextGenerator->Configure( aNode->node_at( nextGenerator->GetName() ) );
            nextGenerator->SetRNG( &fRNG );
            nextGenerator = nextGenerator->GetNextGenerator();
        }

        LMCINFO( lmclog, "Generator toolbox configuration complete" );

        return;
    }

    const Generator* GeneratorToolbox::GetFirstGenerator() const
    {
        return fFirstGenerator;
    }


} /* namespace locust */

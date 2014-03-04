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

    bool GeneratorToolbox::Configure( const ParamNode* aNode )
    {
        if( aNode == NULL ) return false;

        LMCINFO( lmclog, "Creating generators" );

        Factory< Generator >* genFactory = Factory< Generator >::GetInstance();

        const ParamArray* generatorList = aNode->ArrayAt( "generators" );
        if( generatorList == NULL )
        {
            LMCERROR( lmclog, "No generator list was found" );
            return false;
        }

        Generator* lastGenerator = NULL;
        for( ParamArray::const_iterator it = generatorList->Begin(); it != generatorList->End(); ++it )
        {
            if( ! (*it)->IsValue() )
            {
                LMCERROR( lmclog, "Non-value-type array element found in generator-list" );
                continue;
            }

            Generator* newGenerator = genFactory->Create( (*it)->AsValue().Get() );
            if( newGenerator == NULL )
            {
                LMCERROR( lmclog, "Unrecognized generator name: " << (*it)->AsValue().Get() );
                continue;
            }

            if( lastGenerator == NULL )
            {
                fFirstGenerator = newGenerator;
                lastGenerator = newGenerator;
                LMCDEBUG( lmclog, "First generator is <" << fFirstGenerator->GetName() << ">" );
            }
            else
            {
                lastGenerator->SetNextGenerator( newGenerator );
                LMCDEBUG( lmclog, "Adding generator <" << lastGenerator->GetName() << ">" );
            }
        }

        LMCINFO( lmclog, "Configuring generators" );

        Generator* nextGenerator = fFirstGenerator;
        while( nextGenerator != NULL )
        {
            LMCINFO( lmclog, "Configuring generator <" << nextGenerator->GetName() << ">" );
            nextGenerator->Configure( aNode->NodeAt( nextGenerator->GetName() ) );
            nextGenerator->SetRNG( &fRNG );
            nextGenerator = nextGenerator->GetNextGenerator();
        }

        LMCINFO( lmclog, "Generator toolbox configuration complete" );

        return true;
    }

    const Generator* GeneratorToolbox::GetFirstGenerator() const
    {
        return fFirstGenerator;
    }

    Generator* GeneratorToolbox::GetFirstGenerator()
    {
        return fFirstGenerator;
    }


} /* namespace locust */

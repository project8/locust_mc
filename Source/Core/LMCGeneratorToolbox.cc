/*
 * LMCGeneratorToolbox.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#include "LMCGeneratorToolbox.hh"

#include "LMCGenerator.hh"
#include "logger.hh"
#include "LMCParam.hh"

#include "factory.hh"

namespace locust
{
    LOGGER( lmclog, "GeneratorToolbox" );

    GeneratorToolbox::GeneratorToolbox() :
            fFirstGenerator( NULL )
    {
    }

    GeneratorToolbox::~GeneratorToolbox()
    {
        Generator* nextGenerator = fFirstGenerator;
        while( nextGenerator != NULL )
        {
            Generator* thisGenerator = nextGenerator;
            DEBUG( lmclog, "Cleaning up " << thisGenerator->GetName() );
            nextGenerator = thisGenerator->GetNextGenerator();
            delete thisGenerator;
        }
    }

    bool GeneratorToolbox::Configure( const ParamNode* aNode )
    {
        if( aNode == NULL ) return false;

        INFO( lmclog, "Creating generators" );

        scarab::factory< Generator >* genFactory = scarab::factory< Generator >::get_instance();

        const ParamArray* generatorList = aNode->ArrayAt( "generators" );
        if( generatorList == NULL )
        {
            ERROR( lmclog, "No generator list was found" );
            return false;
        } 

        Generator* lastGenerator = NULL;
        for( ParamArray::const_iterator it = generatorList->Begin(); it != generatorList->End(); ++it )
        {
            if( ! (*it)->IsValue() )
            {
                ERROR( lmclog, "Non-value-type array element found in generator-list" );
                continue;
            }
//            else
//            {
//                DEBUG( lmclog, "Reading in reasonable generator: " << (*it)->AsValue().Get() );
//            }

            Generator* newGenerator = genFactory->create( (*it)->AsValue().Get() );
//            DEBUG( lmclog, "And the new generator is ... " << newGenerator->GetName());
            if( newGenerator == NULL )
            {
                ERROR( lmclog, "Unrecognized generator name: " << (*it)->AsValue().Get() );
                continue;
            }

            if( lastGenerator == NULL )
            {
                fFirstGenerator = newGenerator;
                lastGenerator = newGenerator;
                DEBUG( lmclog, "First generator is <" << fFirstGenerator->GetName() << ">" );
            }
            else
            {
//                DEBUG( lmclog, "About to set lastGenerator to ... " << newGenerator->GetName());
                lastGenerator->SetNextGenerator( newGenerator );
                DEBUG( lmclog, "Adding generator <" << lastGenerator->GetNextGenerator()->GetName() << ">");  
                DEBUG( lmclog, "Meanwhile the previous one was <" << lastGenerator->GetName() << ">" );
                lastGenerator = lastGenerator->GetNextGenerator();  // pls addition to advance list pointer.
            }
        }

        INFO( lmclog, "Configuring generators" );

        Generator* nextGenerator = fFirstGenerator;
        while( nextGenerator != NULL )
        {
            INFO( lmclog, "Configuring generator <" << nextGenerator->GetName() << ">" );
            nextGenerator->Configure( aNode->NodeAt( nextGenerator->GetName() ) );
            nextGenerator = nextGenerator->GetNextGenerator();
        }

        INFO( lmclog, "Generator toolbox configuration complete" );

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

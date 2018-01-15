/*
 * LMCGeneratorToolbox.cc
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#include "LMCGeneratorToolbox.hh"

#include "LMCGenerator.hh"
#include "logger.hh"
#include "param.hh"

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
            LDEBUG( lmclog, "Cleaning up " << thisGenerator->GetName() );
            nextGenerator = thisGenerator->GetNextGenerator();
            delete thisGenerator;
        }
    }

    bool GeneratorToolbox::Configure( const scarab::param_node* aNode )
    {
        if( aNode == NULL ) return false;

        LINFO( lmclog, "Creating generators" );

        scarab::factory< Generator >* genFactory = scarab::factory< Generator >::get_instance();

        const scarab::param_array* generatorList = aNode->array_at( "generators" );
        if( generatorList == NULL )
        {
            LERROR( lmclog, "No generator list was found" );
            return false;
        } 

        Generator* lastGenerator = NULL;
        for( scarab::param_array::const_iterator it = generatorList->begin(); it != generatorList->end(); ++it )
        {
            if( ! (*it)->is_value() )
            {
                LERROR( lmclog, "Non-value-type array element found in generator-list" );
                continue;
            }
//            else
//            {
//                LDEBUG( lmclog, "Reading in reasonable generator: " << (*it)->AsValue().Get() );
//            }

            Generator* newGenerator = genFactory->create( (*it)->as_value().as_string() );
//            LDEBUG( lmclog, "And the new generator is ... " << newGenerator->GetName());
            if( newGenerator == NULL )
            {
                LERROR( lmclog, "Unrecognized generator name: " << (*it)->as_value().as_string() );
                continue;
            }

            if( lastGenerator == NULL )
            {
                fFirstGenerator = newGenerator;
                lastGenerator = newGenerator;
                LDEBUG( lmclog, "First generator is <" << fFirstGenerator->GetName() << ">" );
            }
            else
            {
//                LDEBUG( lmclog, "About to set lastGenerator to ... " << newGenerator->GetName());
                lastGenerator->SetNextGenerator( newGenerator );
                LDEBUG( lmclog, "Adding generator <" << lastGenerator->GetNextGenerator()->GetName() << ">");  
                LDEBUG( lmclog, "Meanwhile the previous one was <" << lastGenerator->GetName() << ">" );
                lastGenerator = lastGenerator->GetNextGenerator();  // pls addition to advance list pointer.
            }
        }

        LINFO( lmclog, "Configuring generators" );

        Generator* nextGenerator = fFirstGenerator;
        while( nextGenerator != NULL )
        {
            LINFO( lmclog, "Configuring generator <" << nextGenerator->GetName() << ">" );
            nextGenerator->Configure( aNode->node_at( nextGenerator->GetName() ) );
            nextGenerator = nextGenerator->GetNextGenerator();
        }

        LINFO( lmclog, "Generator toolbox configuration complete" );

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

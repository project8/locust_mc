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

    bool GeneratorToolbox::Configure( const scarab::param_node& aNode )
    {

        LINFO( lmclog, "Creating generators" );

        scarab::factory< Generator >* genFactory = scarab::factory< Generator >::get_instance();

        // TODO: this line will throw an exception if "generators" is not present or it's not an array
        // TODO: this should either check that those are the case and return false if not, or
        // TODO: catch the exception and then return false
        const scarab::param_array& generatorList = aNode["generators"].as_array();




        Generator* lastGenerator = nullptr;

        for( scarab::param_array::const_iterator it = generatorList.begin(); it != generatorList.end(); ++it )
        {
            if( ! it->is_value() )
            {
                LERROR( lmclog, "Non-value-type array element found in generator-list" );
                // TODO: this indicates a problem in the config and should result in locust exiting
                continue;
            }
//            else
//            {
//                LDEBUG( lmclog, "Reading in reasonable generator: " << (*it)->AsValue().Get() );
//            }

            Generator* newGenerator = genFactory->create( (*it)().as_string() );
//            LDEBUG( lmclog, "And the new generator is ... " << newGenerator->GetName());
            if( newGenerator == nullptr )
            {
                LERROR( lmclog, "Unrecognized generator name: " << (*it)().as_string() );
                // TODO: this should also be a fatal error
                continue;
            }

            if( lastGenerator == nullptr )
            {
                fFirstGenerator = newGenerator;
                lastGenerator = newGenerator;
                LINFO( lmclog, "First generator is <" << fFirstGenerator->GetName() << ">" );
            }
            else
            {
//                LDEBUG( lmclog, "About to set lastGenerator to ... " << newGenerator->GetName());
                lastGenerator->SetNextGenerator( newGenerator );
                LINFO( lmclog, "Adding generator <" << lastGenerator->GetNextGenerator()->GetName() << ">");
                LDEBUG( lmclog, "Meanwhile the previous one was <" << lastGenerator->GetName() << ">" );
                lastGenerator = lastGenerator->GetNextGenerator();
            }
        }

        LINFO( lmclog, "Configuring generators" );

        Generator* nextGenerator = fFirstGenerator;
        while( nextGenerator != nullptr )
        {
            LINFO( lmclog, "Configuring generator <" << nextGenerator->GetName() << ">" );
            if( ! aNode.has( nextGenerator->GetName() ) )
            {
                LDEBUG( lmclog, "No configuration information present" );
                continue;
            }
            nextGenerator->Configure( aNode[ nextGenerator->GetName() ].as_node() );
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

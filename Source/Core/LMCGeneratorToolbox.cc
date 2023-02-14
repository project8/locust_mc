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
        // DONE: Since both Scarab v1 and v2 are still being supported, we are returning false
        // DONE: if v1 fails, then we are trying v2, and then throwing an exception if it still
        // DONE: fails.  In the latter case, Scarab presently throws the exception before
        // DONE: Locust does it.
        // TODO: Probably upgrade to Scarab v3.

        if (!(aNode.has("generators")&&aNode["generators"].is_array()))
        {
            LPROG( lmclog, "Trying to parse parameters." );
            // Return false if the parsing failed, for example if trying v1 parsing
            // on a v2 config, or if there is no config, or some other problem.
            return false;
        }

        const scarab::param_array& generatorList = aNode["generators"].as_array();
        if (!(aNode.has("generators")&&aNode["generators"].is_array()))
        {
        	// This exception is not presently thrown, due to Scarab throwing
        	// an exception first.
        	throw std::runtime_error("Parsing either v1 and v2 has not worked.");
            return false;
        }
        else
        {
        	// Everything is working:
        	LPROG( lmclog, "Parsing parameters now.");
        }


        Generator* lastGenerator = nullptr;

        for( scarab::param_array::const_iterator it = generatorList.begin(); it != generatorList.end(); ++it )
        {
            if( ! it->is_value() )
            {
                LERROR( lmclog, "Non-value-type array element found in generator-list" );
                exit(-1);
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
                exit(-1);
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
            	nextGenerator = nextGenerator->GetNextGenerator();
                continue;
            }
            else
            {
            	nextGenerator->Configure( aNode[ nextGenerator->GetName() ].as_node() );
            	nextGenerator = nextGenerator->GetNextGenerator();
            }

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

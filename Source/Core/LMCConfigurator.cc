/*
 * LMCConfigurator.cc
 *
 *  Created on: Nov 5, 2013
 *      Author: nsoblath
 */

#include "LMCConfigurator.hh"

#include "LMCParser.hh"
#include "Logger.hh"

using std::string;

namespace locust
{
    LOGGER( mtlog, "Configurator" );

    Configurator::Configurator( int an_argc, char** an_argv, ParamNode* a_default ) :
            f_master_config( new ParamNode() ),
            f_param_buffer( NULL ),
            f_string_buffer()
    {
        Parser t_parser( an_argc, an_argv );
        //std::cout << "options parsed" << std::endl;
        //cout << t_parser );

        // first configuration: defaults
        if ( a_default != NULL )
        {
            f_master_config->Merge(*a_default);
        }

        //std::cout << "first configuration complete" << std::endl;
        //cout << f_master_config );
        //cout << t_parser );

        string t_name_config("config");
        string t_name_json("json");

        // second configuration: config file
        if( t_parser.Has( t_name_config ) )
        {
            string t_config_filename = t_parser.ValueAt( t_name_config )->Get();
            if( ! t_config_filename.empty() )
            {
                ParamNode* t_config_from_file = ParamInputJSON::ReadFile( t_config_filename );
                if( t_config_from_file == NULL )
                {
                    throw Exception() << "error parsing config file";
                }
                f_master_config->Merge( *t_config_from_file );
                delete t_config_from_file;
            }
        }

        //std::cout << "second configuration complete" << std::endl;
        //cout << f_master_config );
        //cout << t_parser );

        // third configuration: command line json
        if( t_parser.Has( t_name_json ) )
        {
            string t_config_json = t_parser.ValueAt( t_name_json )->Get();
            if( ! t_config_json.empty() )
            {
                ParamNode* t_config_from_json = ParamInputJSON::ReadString( t_config_json );
                f_master_config->Merge( *t_config_from_json );
                delete t_config_from_json;
            }
        }

        //std::cout << "third configuration complete" << std::endl;
        //cout << f_master_config );
        //cout << t_parser );

        // fourth configuration: command line arguments
        t_parser.Erase( t_name_config );
        t_parser.Erase( t_name_json );

        //std::cout << "removed config and json from parsed options" << std::endl;
        //cout << t_parser );
        f_master_config->Merge( t_parser );

        //std::cout << "fourth configuration complete" << std::endl;
        INFO( mtlog, "final configuration:\n" << *f_master_config );
    }

    Configurator::~Configurator()
    {
        delete f_master_config;
    }

    ParamNode* Configurator::Config()
    {
        return f_master_config;
    }

    const ParamNode* Configurator::Config() const
    {
        return f_master_config;
    }

} /* namespace locust */

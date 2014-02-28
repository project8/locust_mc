/*
 * LMCConfigurator.hh
 *
 *  Created on: Nov 5, 2013
 *      Author: nsoblath
 */

#ifndef LMCCONFIGURATOR_HH_
#define LMCCONFIGURATOR_HH_

#include "LMCParam.hh"

#include "LMCException.hh"

#include <string>

namespace locust
{

    class Configurator
    {
        public:
            Configurator( int an_argc, char** an_argv, ParamNode* a_default = NULL );
            virtual ~Configurator();

            ParamNode* Config();
            const ParamNode* Config() const;

            template< typename XReturnType >
            XReturnType Get( const std::string& a_name ) const;

            template< typename XReturnType >
            XReturnType Get( const std::string& a_name, XReturnType a_default ) const;

        private:
            ParamNode* f_master_config;

            mutable Param* f_param_buffer;

            std::string f_string_buffer;
    };

    template< typename XReturnType >
    XReturnType Configurator::Get( const std::string& a_name ) const
    {
        f_param_buffer = const_cast< Param* >( f_master_config->At( a_name ) );
        if( f_param_buffer != NULL && f_param_buffer->IsValue() )
        {
            return f_param_buffer->AsValue().Get< XReturnType >();
        }
        throw Exception() << "Configurator does not have a value for <" << a_name << ">";
    }

    template< typename XReturnType >
    XReturnType Configurator::Get( const std::string& a_name, XReturnType a_default ) const
    {
        f_param_buffer = const_cast< Param* >( f_master_config->At( a_name ) );
        if( f_param_buffer != NULL && f_param_buffer->IsValue() )
        {
            return f_param_buffer->AsValue().Get< XReturnType >();
        }
        return a_default;

    }


} /* namespace locust */
#endif /* LMCCONFIGURATOR_HH_ */

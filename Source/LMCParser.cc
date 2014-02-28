#include "LMCParser.hh"

namespace locust
{

    Parser::Parser( int an_argc, char** an_argv ) :
            ParamNode()
    {
        Parse(an_argc, an_argv);
    }

    Parser::~Parser()
    {
    }

    void Parser::Parse( int an_argc, char** an_argv )
    {
        for( int t_index = 1; t_index < an_argc; t_index++ )
        {
            //t_argument.assign( an_argv[ t_index ] );
            std::string t_argument( an_argv[ t_index ] );
            size_t t_val_pos = t_argument.find_first_of( f_separator );
            if( t_val_pos != std::string::npos )
            {
                std::string t_name(t_argument.substr( 0, t_val_pos ));

                ParamValue* new_value = new ParamValue();
                *new_value << t_argument.substr( t_val_pos + 1 );

                //std::cout << "(Parser) adding < " << t_name << "<" << t_type << "> > = <" << new_value.value() << ">" << std::endl;

                this->Replace( t_name, new_value );

                continue;
            }

            throw Exception() << "argument <" << t_argument << "> does not match <name>=<value> pattern";
        }

        return;
    }


}

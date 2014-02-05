#ifndef MT_PARSER_HH_
#define MT_PARSER_HH_

#include "LMCParam.hh"

#include "LMCException.hh"

#include <map>
#include <string>
#include <sstream>

namespace locust
{

    class Parser : public ParamNode
    {
        public:
            Parser( int an_argc, char** an_argv );
            virtual ~Parser();

            void parse( int an_argc, char** an_argv );

        private:
            static const char f_separator = '=';
            static const size_t f_npos = std::string::npos;

    };

}

#endif

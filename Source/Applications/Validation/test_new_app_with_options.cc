/*
 * test_app_with_options.cc
 *
 *  Created on: Oct 28, 2018
 *      Author: N.S. Oblath
 *
 *  Output examples:
 *
    > bin/test_app_with_options -f 5 -s "hello" -t "world"
    2018-10-29 11:15:20 [ PROG] (tid 0x7fffb0945380) i/application.cc(97): Final configuration:

    {
        third-value : world
    }

    2018-10-29 11:15:20 [ PROG] (tid 0x7fffb0945380) i/application.cc(98): Ordered args:

    [
    ]

    2018-10-29 11:15:20 [ PROG] (tid 0x7fffb0945380) _with_options.cc(69): First value: 5
    2018-10-29 11:15:20 [ PROG] (tid 0x7fffb0945380) _with_options.cc(70): Second value: hello
    2018-10-29 11:15:20 [ PROG] (tid 0x7fffb0945380) _with_options.cc(71): Third value: world
 *
 */

#include "application.hh"

#include "logger.hh"

using namespace scarab;

LOGGER( testlog, "test_app_with_options" );

class test_app : public main_app
{
    public:
        test_app() :
            main_app(),
            f_first_value(),
            f_second_value()
        {
            // Add an option that sets the value directly
            add_option("-f,--first-value", f_first_value, "Set the first value" );

            // Add an option that sets the value via callback
            add_option( "-s,--second-value", [this](std::vector< std::string > args) { f_second_value = args[0]; return true; }, "Set the second value" );

            // Add an option that sets the value into the configuration via callback
            add_option( "-t,--third-value", [this](std::vector< std::string > args) { f_master_config.add( "third-value", args[0] ); return true; }, "Set the third value" );
        }
        virtual ~test_app() {}

        void set_first_value( int a_value ) { f_first_value = a_value; }
        int get_first_value() { return f_first_value; }

        const std::string& second_value() const { return f_second_value; }
        std::string& second_value() { return f_second_value; }

    private:
        int f_first_value;
        std::string f_second_value;
};

int main( int argc, char **argv )
{
    test_app the_main;

    CLI11_PARSE( the_main, argc, argv );

    LPROG( testlog, "First value: " << the_main.get_first_value() );
    LPROG( testlog, "Second value: " << the_main.second_value() );
    LPROG( testlog, "Third value: " << the_main.master_config()["third-value"]() );

    return 0;
}

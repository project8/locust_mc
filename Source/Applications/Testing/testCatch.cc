// testCatch.cc

#define CATCH_CONFIG_RUNNER

#include "application.hh"
#include "logger.hh"
#include "catch.hpp"
#include "LMCTestParameterHandler.hh"


// Define the static Singleton pointer
TestParameterHandler* TestParameterHandler::inst_ = NULL;

TestParameterHandler* TestParameterHandler::getInstance()
{
    if (inst_ == NULL)
    {
    	inst_ = new TestParameterHandler();
    }
   return(inst_);
}



int main(int argc, char *argv[])
{
	TestParameterHandler* p1 = TestParameterHandler::getInstance();
	p1->SetArgs(argc, argv);

	Catch::Session session; // There must be exactly one instance

	return session.run();
}


int Factorial( int number ) {
   return number <= 1 ? number : Factorial( number - 1 ) * number;  // fail
}


TEST_CASE( "Factorials of 1 and higher are computed (pass)", "[single-file]" ) {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}


// Compile & run:
// - g++ -std=c++11 -Wall -I$(CATCH_SINGLE_INCLUDE) -o 010-TestCase 010-TestCase.cpp && 010-TestCase --success
// - cl -EHsc -I%CATCH_SINGLE_INCLUDE% 010-TestCase.cpp && 010-TestCase --success


// testCatch.cc

// Let Catch provide main():
//#define CATCH_CONFIG_MAIN

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
// return number <= 1 ? 1      : Factorial( number - 1 ) * number;  // pass
}

// Mock failed test:
//TEST_CASE( "Factorial of 0 is 1 (fail)", "[single-file]" ) {
//    REQUIRE( Factorial(0) == 1 );
//}

TEST_CASE( "Factorials of 1 and higher are computed (pass)", "[single-file]" ) {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}


// Compile & run:
// - g++ -std=c++11 -Wall -I$(CATCH_SINGLE_INCLUDE) -o 010-TestCase 010-TestCase.cpp && 010-TestCase --success
// - cl -EHsc -I%CATCH_SINGLE_INCLUDE% 010-TestCase.cpp && 010-TestCase --success

// Expected compact output (all assertions):
//
// prompt> 010-TestCase --reporter compact --success
// 010-TestCase.cpp:14: failed: Factorial(0) == 1 for: 0 == 1
// 010-TestCase.cpp:18: passed: Factorial(1) == 1 for: 1 == 1
// 010-TestCase.cpp:19: passed: Factorial(2) == 2 for: 2 == 2
// 010-TestCase.cpp:20: passed: Factorial(3) == 6 for: 6 == 6
// 010-TestCase.cpp:21: passed: Factorial(10) == 3628800 for: 3628800 (0x375f00) == 3628800 (0x375f00)
// Failed 1 test case, failed 1 assertion.

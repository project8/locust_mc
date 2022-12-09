/*
 * LMCTestParameterHandler.hh
 *
 *  Created on: Dec 7, 2022
 *      Author: pslocum
 */

#ifndef LMCTESTPARAMETERHANDLER_HH_
#define LMCTESTPARAMETERHANDLER_HH_

#include <iostream>
using namespace std;

 /*!
 @class TestParameterHandler
 @author P. Slocum
 @brief Class to configure Catch2 test cases.
 @details
 Available configuration options:
 No input parameters
 */

class TestParameterHandler
{

    public:
	    static TestParameterHandler* getInstance();
	    void SetParameters(double aValue) {testVar = aValue;};
	    double GetParameters() {return(testVar);};
	    void SetArgs(int anArgc, char** anArgv) {argc = anArgc; argv = anArgv;};
	    int GetArgc() {return(argc);};
	    char** GetArgv() {return(argv);};



    private:
       static TestParameterHandler* inst_;   // The one, single instance
       TestParameterHandler() : testVar(0), argc(0), argv(0) {} // private constructor
       TestParameterHandler(const TestParameterHandler&);
       TestParameterHandler& operator=(const TestParameterHandler&);

       double testVar;
       int argc;
       char** argv;
};




#endif /* LMCTESTPARAMETERHANDLER_HH_ */

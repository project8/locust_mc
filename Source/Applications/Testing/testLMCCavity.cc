/*
 *  testLMCCavity.cc
 *
 *  Created on: Jan 24, 2023
 *      Author: P. L. Slocum
 *
 *  This unit test calculates an unnormalized response of a Green's function
 *  corresponding to a damped resonant harmonic oscillator.  By driving the
 *  oscillator with a stepped range of frequencies around the resonance, the
 *  Q of the resonance can be checked as Q = f_res / Delta_f, where f_res
 *  inferred resonant frequency and Delta_f is the inferred FWHM of the
 *  resonance.  This may be helpful in optimizing parameterized cavity settings
 *  used in the larger simulation.
 *
 *  Tunable parameters include:
 *
 *  "dho-cavity-resonance" [1.067e9 Hz] Resonant frequency of the damped harmonic oscillator (dho).
 *  "dho-cavity-Q" [1000] Q of the dho.
 *  "dho-time-resolution" [1.e-8 s] Bin width of the FIR filter derived from the Green's function.
 *  "dho-threshold-factor" [0.01] Cut-off threshold for FIR ring-down amplitude, relative to max.
 *
 *  It runs like this:
 *  /path/to/testLMCCavity
 *
 *  And should show this final output:
 *  2023-01-24 10:53:25 [ PROG] (tid 140342977077888) testLMCCavity.cc(224): Estimated Q is 1000
 *  2023-01-24 10:53:25 [ PROG] (tid 140342977077888) testLMCCavity.cc(225): Expected Q is 1000
 *
 *  We can try changing a parameter like this:
 *  /path/to/testLMCCavity -m 0.1
 *
 *  And then check the final output:
 *  2023-01-24 11:00:27 [ PROG] (tid 140091322118784) testLMCCavity.cc(224): Estimated Q is 714.286
 *  2023-01-24 11:00:27 [ PROG] (tid 140091322118784) testLMCCavity.cc(225): Expected Q is 1000
 *
 *  The above shows how the inferred Q of the oscillator responds to e.g. the "dho-threshold-factor"
 *  parameter.  Similar tests with the other parameters can be used to evaluate e.g. fidelity of the
 *  resonance, along with other metrics such as speed.  To display the list of available parameters,
 *  /path/to/testLMCCavity -h
 *
 */

#include "LMCSignal.hh"
#include "LMCTFFileHandler.hh"
#include "LMCFIRFileHandler.hh"
#include "LMCAnalyticResponseFunction.hh"
#include "LMCDampedHarmonicOscillator.hh"
#include "application.hh"
#include "logger.hh"
#include "LMCConst.hh"
#include <fftw3.h>
#include <math.h>
#include "catch.hpp"
#include "LMCTestParameterHandler.hh"
#include "LMCCavityUtility.hh"

using namespace scarab;
using namespace locust;

LOGGER( testlog, "testLMCCavity" );

class testCavity_app : public main_app
{
    public:
        testCavity_app() :
        main_app(),
        fDHOTimeResolution(1.e-8),
        fDHOThresholdFactor(0.01),
        fCavityFrequency(1.067e9),
        fCavityQ(1000.),
        fExpandSweep(1.0),
        fUnitTestOutputFile(false),
        fOutputPath( TOSTRING(PB_OUTPUT_DIR) )
        {
            add_option("-r,--dho-time-resolution", fDHOTimeResolution, "[1.e-8] Time resolution used in Green's function (s).");
            add_option("-m,--dho-threshold-factor", fDHOThresholdFactor, "[0.01] Minimum fractional threshold of Green's function used to calculate FIR.");
            add_option("-f,--dho-cavity-frequency", fCavityFrequency, "[1.067e9] Cavity resonant frequency (Hz).");
            add_option("-g,--dho-cavity-Q", fCavityQ, "[1000] Cavity Q.");
            add_option("-x,--expand-sweep", fExpandSweep, "[1.0] Factor by which to expand range of frequency sweep.");
            add_option("-w, --write-output", fUnitTestOutputFile, "[0==false] Write histo to Root file.");
            add_option("-o, --output-path", fOutputPath, "[PB_OUTPUT_DIR]");
        }

        virtual ~testCavity_app() {}

        double GetDHOTimeResolution()
        {
        	return fDHOTimeResolution;
        }
        double GetDHOThresholdFactor()
        {
        	return fDHOThresholdFactor;
        }
        double GetCavityFrequency()
        {
        	return fCavityFrequency;
        }
        double GetCavityQ()
        {
        	return fCavityQ;
        }
        double GetExpandSweep()
        {
        	return fExpandSweep;
        }
        bool UnitTestOutputFile()
        {
        	return fUnitTestOutputFile;
        }
        std::string GetOutputPath()
        {
        	return fOutputPath;
        }


    private:
        double fDHOTimeResolution;
        double fDHOThresholdFactor;
        double fCavityFrequency;
        double fCavityQ;
        double fExpandSweep;
        bool fUnitTestOutputFile;
        std::string fOutputPath;
};


bool parseCavity(testCavity_app& the_main)
{
	TestParameterHandler* p1 = TestParameterHandler::getInstance();
	CLI11_PARSE( the_main, p1->GetArgc(), p1->GetArgv() );
	return true;
}


TEST_CASE( "testLMCCavity with default parameter values (pass)", "[single-file]" )
{
	testCavity_app the_main;
	if (!parseCavity(the_main)) exit(-1);

	CavityUtility aCavityUtility;

	aCavityUtility.SetOutputPath(the_main.GetOutputPath());
	aCavityUtility.SetExpandFactor(the_main.GetExpandSweep());
	aCavityUtility.SetOutputFile(the_main.UnitTestOutputFile());
	// TO-DO:  Remove hard-wired TE011 mode in the next line, and make it more flexible for other modes.
	bool checkCavityQ = aCavityUtility.CheckCavityQ( 1, 0, 1, 1, the_main.GetDHOTimeResolution(), the_main.GetDHOThresholdFactor(), the_main.GetCavityFrequency(), the_main.GetCavityQ() );

	REQUIRE( checkCavityQ );

}



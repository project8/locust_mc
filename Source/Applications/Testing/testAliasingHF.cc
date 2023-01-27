/*
 *  testAliasingHF.cc
 *
 *  Created on: Jan 26, 2023
 *      Author: P. L. Slocum
 *
 *  This unit test calculates HF aliasing.
 *
 *  Tunable parameters include:
 *
 *  "lo-frequency": Local oscillator frequency.
 *
 *
 *  It runs like this:
 *  /path/to/testAliasingHF
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
#include "application.hh"
#include "logger.hh"
#include "catch.hpp"
#include "LMCTestParameterHandler.hh"

using namespace scarab;
using namespace locust;

LOGGER( testlog, "testAliasingHF" );

class testAliasingHF_app : public main_app
{
    public:
        testAliasingHF_app() :
            main_app(),
			fLOFrequency(25.8802e9),
			fRFFrequency(25.9002e9),
			fAcquisitionRate(205.),
			fDecimationRate(10.)
        {
            add_option("-l,--lo-frequency", fLOFrequency, "[25.8802e9] Frequency of local oscillator (Hz).");
            add_option("-r,--rf-frequency", fRFFrequency, "[25.9002e9] RF frequency (Hz).");
            add_option("-a,--acquisition-rate", fAcquisitionRate, "[205.] Data acquisition rate (MHz).");
            add_option("-d,--decimation-rate", fDecimationRate, "[10.] Decimation factor (hard-wired to 10 in Locust.");
        }

        virtual ~testAliasingHF_app() {}

        double GetLocalOscillator()
        {
        	return fLOFrequency;
        }
        double GetRFFrequency()
        {
        	return fRFFrequency;
        }
        double GetAcquisitionRate()
        {
        	return fAcquisitionRate;
        }
        double GetDecimationRate()
        {
        	return fDecimationRate;
        }

    private:
        double fLOFrequency;
        double fRFFrequency;
        double fAcquisitionRate;
        double fDecimationRate;
};


class testAliasingHF
{
public:

	testAliasingHF():
        fLOFrequency(0.),
		fRFFrequency(0.),
        fAcquisitionRate(0.),
		fDecimationRate(0.),
	    param_0( new param_node() )
    {
    }

	void AddParam(std::string aString, double aValue)
	{
		if ( aValue != 0. )
		{
            param_0->add(aString.c_str(), aValue);
		}
	}

	const scarab::param_node* GetParams()
	{
	    return param_0;
	}

	bool Configure()
	{
		return true;
	}


    scarab::param_node* param_0;


    double fLOFrequency;
    double fRFFrequency;
    double fAcquisitionRate;
    double fDecimationRate;

};


bool parseAliasing(testAliasingHF_app& the_main)
{
	TestParameterHandler* p1 = TestParameterHandler::getInstance();
	CLI11_PARSE( the_main, p1->GetArgc(), p1->GetArgv() );
	return true;
}


TEST_CASE( "testAliasingHF with default parameter values (pass)", "[single-file]" )
{
	testAliasingHF_app the_main;
	if (!parseAliasing(the_main)) exit(-1);

	testAliasingHF aTestAliasingHF;
	aTestAliasingHF.AddParam( "lo-frequency", the_main.GetLocalOscillator() );

	if (!aTestAliasingHF.Configure())
	{
		LWARN(testlog,"testAliasingHF was not configured correctly.");
	    REQUIRE( 0 > 1 );
	}

	int sign = 1; // start with upper sideband.
	double LO = the_main.GetLocalOscillator();
	double fs = 1.e6 * the_main.GetAcquisitionRate() * the_main.GetDecimationRate();
	double RF = the_main.GetRFFrequency();

    LPROG( testlog, "\nSummary:");
	LPROG( testlog, "Presently hard-wired decimation rate is 10." );
	LPROG( testlog, "acquisition-rate is " << the_main.GetAcquisitionRate() );
	LPROG( testlog, "lo-frequency is " << the_main.GetLocalOscillator() );
	LPROG( testlog, "rf-frequency is " << the_main.GetRFFrequency() );


	bool bPass = true;

	for (int j=0; j<2; j++)  // choose upper or lower mixing sideband
	{
	    if (j==1) sign = -1;
	    for (int i=0; i<3; i++)
	    {
	    double freq = i*RF - sign*LO;
	    int nwin = (round)(freq/fs);

	    double alias = sign*(freq - nwin*fs);
	    if (j==0)
	    {
	    	printf("%5d*RF - LO = %10.4g and alias is %2d*(%10.4g - %4d*%g)=%10.4g Hz\n", i, freq, sign, freq, nwin, fs, alias);
	    	if ((i!=1) && (fabs(alias) < the_main.GetAcquisitionRate()*1.e6/2.))
	    	{
	    		LERROR( testlog, "Aliased frequency " << fabs(alias) << " is below Nyquist frequency " << the_main.GetAcquisitionRate()*1.e6/2. );
	    		bPass = false;
	    	}
	    }
	    if (j==1)
	    {
	    	printf("%5d*RF + LO = %10.4g and alias is %2d*(%10.4g - %4d*%g)=%10.4g Hz\n", i, freq, sign, freq, nwin, fs, alias);
	    	if (fabs(alias) < the_main.GetAcquisitionRate()*1.e6/2.)
	    	{
	    		LERROR( testlog, "Aliased frequency " << fabs(alias) << " is below Nyquist frequency " << the_main.GetAcquisitionRate()*1.e6/2. );
	    		bPass = false;
	    	}
	    }

	    } // i
	} // j


    REQUIRE( bPass );
}



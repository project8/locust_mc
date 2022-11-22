// testCatchLocust
//#define CATCH_CONFIG_MAIN


#include "application.hh"
#include "logger.hh"
#include <fftw3.h>
#include <math.h>
#include "catch.hpp"

using namespace scarab;

LOGGER( testlog, "testFFT" );

class test_app : public main_app
{
    public:
        test_app() :
            main_app()
        {
        }

        virtual ~test_app() {}

    private:

};


double GetPower(bool bTime)
{
	fftw_plan plan;
	fftw_complex *data;
	int N0 = 1000;
	/* create plan for forward DFT */
	data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N0 );
	plan = fftw_plan_dft_1d(N0, data, data, FFTW_FORWARD, FFTW_ESTIMATE);

	/* initialize data */
	double pdata = 0.;
	for (int j = 0; j < N0; ++j)
	{
	    // IQ signal:
	    data[j][0] = 1.0 * sin(-0.1*j);
	    data[j][1] = 1.0 * cos(0.1*j);  // (set this to 0. for real signal.)

	    // power in time series:
	    pdata += data[j][0]*data[j][0]+data[j][1]*data[j][1];
	}

	/* compute transform, in-place */
	fftw_execute(plan);

	double normalization = sqrt((double)N0);
	double ptransform = 0.;
	for (int j = 0; j < N0; ++j)
	{
	    // IQ signal:
		data[j][0] /= normalization;
		data[j][1] /= normalization;

		// power in frequency spectrum:
	    ptransform += data[j][0]*data[j][0]+data[j][1]*data[j][1];
	}

    if (bTime)
    {
    	LPROG(testlog, "power of original data time series is: " << pdata);
    	return pdata;
    }
    else
    {
  	    LPROG(testlog, "power of transform is: " << ptransform);
    	return ptransform;
    }

}



TEST_CASE( "Power before FFT = power after FFT (pass)", "[single-file]" )
{
	double threshold = 1.e-4;
    REQUIRE( fabs(GetPower(1)-GetPower(0)) < threshold );
}



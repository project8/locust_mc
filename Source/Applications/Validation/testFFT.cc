/*
 * testFFT.cc
 *
 *  Created on: Jun 6, 2022
 *      Author: P. L. Slocum
 *
 *  Output examples:
 *
 */

#include "application.hh"
#include "logger.hh"
#include <fftw3.h>
#include <math.h>

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

int main( int argc, char **argv )
{

    test_app the_main;
    CLI11_PARSE( the_main, argc, argv );

	fftw_plan plan;
	fftw_complex *data;
	int N0 = 1000;
	/* create plan for forward DFT */
	data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N0 );
	plan = fftw_plan_dft_1d(N0, data, data, FFTW_FORWARD, FFTW_ESTIMATE);

	/* initialize data */
	double pdata = 0;
	for (int j = 0; j < N0; ++j)
	{
	    // IQ signal:
	    data[j][0] = 1.0 * sin(-0.1*j);
	    data[j][1] = 1.0 * cos(0.1*j);  // (set this to 0. for real signal.)

	    // power in time series:
	    pdata += data[j][0]*data[j][0]+data[j][1]*data[j][1];
	}

	LPROG(testlog, "power of original data time series is: " << pdata);

	/* compute transform, in-place */
	fftw_execute(plan);

	double normalization = sqrt((double)N0);
	double ptransform = 0;
	for (int j = 0; j < N0; ++j)
	{
	    // IQ signal:
		data[j][0] /= normalization;
		data[j][1] /= normalization;

		// power in frequency spectrum:
	    ptransform += data[j][0]*data[j][0]+data[j][1]*data[j][1];
	}


	  LPROG(testlog, "power of transform is: " << ptransform);

	  if (fabs(pdata-ptransform) < pdata*1.e-4)
	  {
		  LPROG(testlog, "Hooray, energy is conserved! ");
		  LPROG(testlog, "Press return to continue ... ");
		  getchar();
	  }
	  else
	  {
		  LERROR(testlog, "Something went wrong, energy is not being conserved.");
		  LERROR(testlog, "power of original data time series is: " << pdata);
		  LERROR(testlog, "but power of transform is: " << ptransform);
		  LPROG(testlog, "Press return to continue ... ");
		  getchar();
	  }


	  fftw_destroy_plan(plan);


    return 0;
}

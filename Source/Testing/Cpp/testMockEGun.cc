/*
 * testMockEGun.cc
 *
 *  Created on: Jun 7, 2022
 *      Author: P. L. Slocum
 *
 *  Output examples:
 *
 *  This unit test calculates power detected in a rectangular waveguide in which a
 *  26 GHz electron is undergoing cyclotron motion.  Total power emitted by the electron
 *  is 1.e-15 W.  One-way power detected at one end of the waveguide is reduced by half,
 *  and again by 0.4 due to average dot product with the mode field.
 *
 *  For example, the command:
 *
 *  > bin/testMockEGun
 *
 *  produces the output below:
 *
 *  2022-06-07 17:08:54 [ PROG] (tid 140543355190208) /testMockEGun.cc(63): power of original data time series is: 2e-16
 *  2022-06-07 17:08:54 [ PROG] (tid 140543355190208) /testMockEGun.cc(81): power of transformed data is: 2e-16
 *
 *  This output aligns with our typical e-gun simulation results.
 */


#include "application.hh"
#include "logger.hh"
#include <fftw3.h>
#include <math.h>
#include "catch.hpp"

using namespace scarab;

LOGGER( testlog, "testMockEGun" );

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

double GetPower()
{

	double larmorPower = 1.e-15;
	double halfPower = 0.5; // Fraction that propagates in one direction in the waveguide.
	double dotProductFactor = 0.4; // Average loss due to dot product between velocity and TE01 in rectangular waveguide.


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
	    data[j][0] = sqrt(larmorPower*halfPower*dotProductFactor) * sin(-0.1*j);
	    data[j][1] = sqrt(larmorPower*halfPower*dotProductFactor) * cos(0.1*j);  // (set this to 0. for real signal.)

	    // power in time series:
	    pdata += data[j][0]*data[j][0]+data[j][1]*data[j][1];
	}

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

	fftw_destroy_plan(plan);

    return ptransform/N0;

}


TEST_CASE( "Larmor power fraction. (pass)", "[single-file]" )
{
	double expectedPower = 2.e-16;
	double threshold = 1.e-4;
    REQUIRE( fabs(GetPower() - 2.e-16) <= threshold*expectedPower );
}




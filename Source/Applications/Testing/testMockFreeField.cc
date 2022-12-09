/*
 * testMockFreeField.cc
 *
 *  Created on: Jun 7, 2022
 *      Author: P. L. Slocum
 *
 *  This unit test calculates total detected power around a point source with total emitted power 1.e-15 W.
 *  As a separate functionality, it also accepts power [Watts] detected at an antenna with an effective
 *  aperture of 1.e-5 m^2, and derives the expected E-field amplitudes at that antenna.
 *
 *  For example, for a detected power of 2.e-19 Watts at the antenna, the command:
 *
 *  bin/testMockFreeField -i 2.e-19
 *
 *  produces the output below, showing that the E-field amplitude expected at the antenna is 2.74e-6 V/m.
 *
 *  2022-06-07 16:57:45 [ PROG] (tid 140283443030976) MockFreeField.cc(90): power of original data is: 1e-15 Watts
 *  2022-06-07 16:57:45 [ PROG] (tid 140283443030976) MockFreeField.cc(108): power of transformed data is: 1e-15 Watts
 *  2022-06-07 16:57:45 [ PROG] (tid 140283443030976) MockFreeField.cc(112): Hooray, energy is conserved!
 *  2022-06-07 16:57:45 [ PROG] (tid 140283443030976) MockFreeField.cc(114): E-field amplitude should be: 2.74492e-06
 *
 *
 */

#include "application.hh"
#include "logger.hh"
#include <fftw3.h>
#include <math.h>
#include "catch.hpp"
#include "LMCConst.hh"
#include "LMCTestParameterHandler.hh"


using namespace scarab;

LOGGER( testlog, "testMockFreeField" );

class test_app : public main_app
{
    public:
        test_app() :
            main_app(),
			fTestParameter(0.)
        {
            add_option("-t,--test-parameter", fTestParameter, "Set a test parameter." );
        }

        virtual ~test_app() {}

        double GetTestParameter()
        {
            return fTestParameter;
        }

    private:
        double fTestParameter;
};


class PowerHandler
{
    public:
    PowerHandler() :
        fLarmorPower( 1.e-15 ), // Watts
        fIncidentPower( 0. ), // Watts
        fRadius( 0.1 ), // m
        fEffectiveAperture( 1.e-5 ) // m^2.
    {
    }

        void SetLarmorPower( double aPower ) { fLarmorPower = aPower; }
        double GetLarmorPower() { return fLarmorPower; }
        void SetIncidentPower( double aPower ) { fIncidentPower = aPower; }
        double GetIncidentPower() { return fIncidentPower; }
        void SetRadius( double aRadius ) { fRadius = aRadius; }
        double GetRadius() { return fRadius; }
        void SetEffectiveAperture( double anArea ) { fEffectiveAperture = anArea; }
        double GetEffectiveAperture() { return fEffectiveAperture; }

    private:

        double fLarmorPower; // Watts
        double fIncidentPower; // Watts
        double fRadius; // m
        double fEffectiveAperture; // m^2.

};

int parseFreeField(test_app& the_main)
{
	TestParameterHandler* p1 = TestParameterHandler::getInstance();
    CLI11_PARSE( the_main, p1->GetArgc(), p1->GetArgv() );
	return 0;
}



TEST_CASE( "Mock free space Larmor power. (pass)", "[single-file]" )
{
	test_app the_main;
	parseFreeField(the_main);
	PowerHandler aPowerHandler;

	double radius = aPowerHandler.GetRadius(); // meters
	double PoyntingVector = aPowerHandler.GetLarmorPower() / (4.*locust::LMCConst::Pi()*radius*radius);
	double isotropicIncidentPower = PoyntingVector * aPowerHandler.GetEffectiveAperture();
	if (aPowerHandler.GetIncidentPower() == 0.)
		{
			double incidentPower = isotropicIncidentPower;
			aPowerHandler.SetIncidentPower(incidentPower);
		}

	fftw_plan plan;
	fftw_complex *data;
	int N0 = 1000;
//	 create plan for forward DFT
	data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N0 );
	plan = fftw_plan_dft_1d(N0, data, data, FFTW_FORWARD, FFTW_ESTIMATE);

//	 initialize data
	double pdata = 0;
	for (int j = 0; j < N0; ++j)
	{
	    // IQ signal:
	    data[j][0] = sqrt(isotropicIncidentPower) * sin(-0.1*j);
	    data[j][1] = sqrt(isotropicIncidentPower) * cos(0.1*j);  // (set this to 0. for real signal.)

	    // power in time series:
	    pdata += data[j][0]*data[j][0]+data[j][1]*data[j][1];
	}

	LPROG(testlog, "power of original data is: " << pdata * 4. * locust::LMCConst::Pi() * radius * radius / aPowerHandler.GetEffectiveAperture() / N0 << " Watts");

//	 compute transform, in-place
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


	  LPROG(testlog, "power of transformed data is: " << ptransform  * 4. * locust::LMCConst::Pi() * radius * radius / aPowerHandler.GetEffectiveAperture() / N0 << " Watts");

	  if (fabs(pdata-ptransform) < 1.e-4*pdata)
	  {
		  LPROG(testlog, "Hooray, energy is conserved! ");
		  // E-amplitude = sqrt( S/(c*epsilon0) )
		  LPROG(testlog, "E-field amplitude should be: " << sqrt(aPowerHandler.GetIncidentPower() / aPowerHandler.GetEffectiveAperture()/locust::LMCConst::C()/locust::LMCConst::EpsNull()));
	  }
	  else
	  {
		  LERROR(testlog, "Something went wrong, energy is not being conserved.");
		  LERROR(testlog, "power of original data time series is: " << pdata);
		  LERROR(testlog, "but power of transform is: " << ptransform);
	  }


	  fftw_destroy_plan(plan);


	double threshold = 1.e-4;
    REQUIRE( fabs( pdata-ptransform) <= threshold*pdata );


}

#include "LMCSignal.hh"
#include "application.hh"
#include "logger.hh"
#include "LMCConst.hh"
#include <fftw3.h>
#include <math.h>
#include "catch.hpp"
#include "LMCTestParameterHandler.hh"


using namespace scarab;
using namespace locust;

LOGGER( testlog, "testLMCTestSignal" );

class test_app : public main_app
{
    public:
        test_app() :
            main_app(),
			fTestParameter(0.)
        {
            add_option("-t,--test-parameter", fTestParameter, "Set a test parameter.");
        }

        virtual ~test_app() {}

        double GetTestParameter()
        {
            return fTestParameter;
        }


    private:
        double fTestParameter;
};


class testLMCTestSignal
{
public:

	testLMCTestSignal():
        fRF_frequency(20.1e9),
        fLO_frequency(20.05e9),
        fAmplitude(5.e-8),
        fAcquisitionRate(200.)
    {
    }


	double GetAmplitude()
	{
		return fAmplitude;
	}

	bool PopulateSignal(Signal* aSignal)
	{

        double LO_phase = 0.;
        double voltage_phase = 0.;
        int nChannels = 1;

        for (unsigned ch = 0; ch < nChannels; ++ch)
        {
            for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
            {

                LO_phase = 2.*LMCConst::Pi()*fLO_frequency*(double)index/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
                voltage_phase = 2.*LMCConst::Pi()*fRF_frequency*(double)index/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);

                // keep only lower sideband RF-LO
                aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0] += sqrt(50.)*fAmplitude*cos(voltage_phase-LO_phase);
                aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1] += sqrt(50.)*fAmplitude*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);

            }
        }
        return true;
	}

    double GetPower(Signal* aSignal, int N0, bool bTime)
    {


	    fftw_plan plan;
	    fftw_complex *data;

	    /* create plan for forward DFT */
	    data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N0 );
	    plan = fftw_plan_dft_1d(N0, data, data, FFTW_FORWARD, FFTW_ESTIMATE);

	    /* initialize data */
	    double pdata = 0.;
	    for (int j = 0; j < N0; ++j)
	    {
	        // IQ signal:
	        data[j][0] = aSignal->LongSignalTimeComplex()[j][0];
	        data[j][1] = aSignal->LongSignalTimeComplex()[j][1];  // (set this to 0. for real signal.)

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
    	    LPROG(testlog, "fft power of original data time series is: " << pdata);
    	    return pdata;
        }
        else
        {
  	        LPROG(testlog, "power of transform is: " << ptransform);
    	    return ptransform;
        }

    }


    private:

    double fRF_frequency;
    double fLO_frequency;
    double fAcquisitionRate;
    double fAmplitude;

};

int parseTestSignal()
{
	test_app the_main;
	TestParameterHandler* p1 = TestParameterHandler::getInstance();
    CLI11_PARSE( the_main, p1->GetArgc(), p1->GetArgv() );
	return 0;
}





TEST_CASE( "LMCTestSignal with default parameter values (pass)", "[single-file]" )
{
	parseTestSignal();
	testLMCTestSignal aTestLMCTestSignal;
    Signal* aSignal = new Signal();
    int N0 = 1000;
    aSignal->Initialize( N0 , 1 );

   /* initialize data */
    aTestLMCTestSignal.PopulateSignal(aSignal);

	double threshold = 1.e-4;
    REQUIRE( fabs( aTestLMCTestSignal.GetPower(aSignal, N0, 1) - aTestLMCTestSignal.GetPower(aSignal, N0, 0) ) < threshold*aTestLMCTestSignal.GetPower(aSignal, N0, 1) );
    REQUIRE( fabs( aTestLMCTestSignal.GetPower(aSignal, N0, 1)/N0/50. - pow(aTestLMCTestSignal.GetAmplitude(),2.) ) < threshold*aTestLMCTestSignal.GetPower(aSignal, N0, 1)/N0/50.);
}



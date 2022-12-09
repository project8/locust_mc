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

using namespace scarab;
using namespace locust;

LOGGER( testlog, "testLMCCavity" );

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


class testLMCCavity
{
public:

	testLMCCavity():
        fRF_frequency(0.),
        fLO_frequency(0.),
        fAmplitude(0.),
        fAcquisitionRate(0.),
		fFilterRate(0.),
		fPhaseLagLO(0.),
		fTFReceiverHandler( 0 ),
		fAnalyticResponseFunction( 0 )
    {
    }

	const scarab::param_node* GetParams()
	{
	    scarab::param_node* param_0( new param_node() );
	    std::string tDataDir = TOSTRING(PB_DATA_INSTALL_DIR);
	    std::string filepath = tDataDir + "/PatchTFLocust.txt";
	    param_0->add( "tf-receiver-filename", filepath.c_str() );
	    param_0->add( "tf-receiver-bin-width", 0.01e9 );
	    return param_0;
	}

	bool Configure()
	{
		fTFReceiverHandler = new TFReceiverHandler();
		if ( !fTFReceiverHandler->Configure(*GetParams()) )
		{
			LWARN(testlog,"TFReceiverHandler was not configured correctly.");
		    return false;
		}

		fAnalyticResponseFunction = new DampedHarmonicOscillator();
		if ( !fAnalyticResponseFunction->Configure(*GetParams()) )
		{
			LWARN(testlog,"DampedHarmonicOscillator was not configured.");
			return false;
		}
		if ( !fTFReceiverHandler->ConvertAnalyticGFtoFIR(fAnalyticResponseFunction->GetGFarray()) )
		{
			LWARN(testlog,"GF->FIR was not generated.");
			return false;
		}

        fRF_frequency = 0.9e9; // Hz
        fLO_frequency = 0.05e9; // Hz
        fAcquisitionRate = 201.e6; // Hz
        fAmplitude = 5.e-6; // volts
        fPhaseLagLO = 0.;
		return true;
	}



	bool PopulateSignal(Signal* aSignal, int N0, double timeStamp)
	{

        double LO_phase = 0.;
        double voltage_phase = 0.;
        double phaseLagRF = 2.*LMCConst::Pi() * timeStamp/(1./fRF_frequency);
        fPhaseLagLO += 2.*LMCConst::Pi() * fLO_frequency * 1./fAcquisitionRate;

        for( unsigned index = 0; index < N0; ++index )
        {
            LO_phase = fPhaseLagLO;
            voltage_phase = phaseLagRF + 2.*LMCConst::Pi()*fRF_frequency*(double)index/fFilterRate;

            // keep only lower sideband RF-LO
            aSignal->LongSignalTimeComplex()[index][0] = fAmplitude*cos(voltage_phase-LO_phase);
            aSignal->LongSignalTimeComplex()[index][1] = fAmplitude*cos(-LMCConst::Pi()/2. + voltage_phase-LO_phase);
        }
//        printf("fAmplitude is %g, signal[0][0] is %g with phaseLagRF %g and LOPhase %g at time stamp %g\n", fAmplitude, aSignal->LongSignalTimeComplex()[0][0], phaseLagRF, fPhaseLagLO, timeStamp);
//        getchar();
        return true;
	}

	std::deque<double> SignalToDeque(Signal* aSignal)
	{
	    std::deque<double> incidentSignal;
	    for (unsigned i=0; i<fTFReceiverHandler->GetFilterSize(); i++)
	    {
	    	incidentSignal.push_back(aSignal->LongSignalTimeComplex()[i][0]);
	    }
	    return incidentSignal;
	}


    TFReceiverHandler* fTFReceiverHandler;
    AnalyticResponseFunction* fAnalyticResponseFunction;

    double fRF_frequency;
    double fLO_frequency;
    double fAcquisitionRate;
    double fFilterRate;
    double fAmplitude;
    double fPhaseLagLO;



};


int parseCavity(test_app& the_main)
{
	TestParameterHandler* p1 = TestParameterHandler::getInstance();
    CLI11_PARSE( the_main, p1->GetArgc(), p1->GetArgv() );
	return 0;
}


TEST_CASE( "testLMCCavity with default parameter values (pass)", "[single-file]" )
{
	test_app the_main;
	parseCavity(the_main);
	LPROG( testlog, "fTestParameter is " << the_main.GetTestParameter() );

	testLMCCavity aTestLMCCavity;
	if (!aTestLMCCavity.Configure())
	{
		LWARN(testlog,"testLMCCavity was not configured correctly.");
	    REQUIRE( 0 > 1 );
	}

    /* initialize time series */
    Signal* aSignal = new Signal();
    int N0 = aTestLMCCavity.fTFReceiverHandler->GetFilterSize();
    aTestLMCCavity.fFilterRate = (1./aTestLMCCavity.fTFReceiverHandler->GetFilterResolution());
    aSignal->Initialize( N0 , 1 );

    double firGainMax = 0.;

        for (unsigned rfStep=0; rfStep<30; rfStep++) // frequency sweep
        {
            aTestLMCCavity.fRF_frequency = 1.0e9 + 0.005e9*rfStep;

	        double convolutionMag = 0.;
            for (unsigned i=0; i<1000; i++)  // time stamps
            {
                // populate time series and convolve it with the FIR filter
        	    double timeStamp = i/aTestLMCCavity.fAcquisitionRate;
                aTestLMCCavity.PopulateSignal(aSignal, N0, timeStamp);
            	std::pair<double,double> convolutionPair = aTestLMCCavity.fTFReceiverHandler->ConvolveWithComplexFIRFilter(aTestLMCCavity.SignalToDeque(aSignal));
                if (fabs(convolutionPair.first) > convolutionMag)
                {
        	        convolutionMag = convolutionPair.first;
                }
            } // i
            double firGain = convolutionMag;
			LPROG( testlog, "Cavity GF gain at frequency " << aTestLMCCavity.fRF_frequency << " is " << firGain );
        } // rfStep

    delete aSignal;

    REQUIRE( 1 > 0. ); // to-do:  get this defined.
}



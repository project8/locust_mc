#include "LMCSignal.hh"
#include "LMCTFFileHandler.hh"
#include "LMCFIRFileHandler.hh"
#include "application.hh"
#include "logger.hh"
#include "LMCConst.hh"
#include <fftw3.h>
#include <math.h>
#include "catch.hpp"
#include "LMCTestParameterHandler.hh"


using namespace scarab;
using namespace locust;

LOGGER( testlog, "testLMCPlaneWaveFIR" );

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


class testLMCPlaneWaveFIR
{
public:

	testLMCPlaneWaveFIR():
        fRF_frequency(0.),
        fLO_frequency(0.),
        fAmplitude(0.),
        fAcquisitionRate(0.),
		fFilterRate(0.),
		fPhaseLagLO(0.),
		fTFReceiverHandler( 0 )
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
		if ( !fTFReceiverHandler->ReadHFSSFile() )
	    {
			LWARN(testlog,"TF file was not read correctly.");
			return false;
	    }

        fRF_frequency = 25.9e9; // Hz
        fLO_frequency = 25.85e9; // Hz
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
//        printf("signal[0][0] is %g with phaseLagRF %g and LOPhase %g at time stamp %g\n", aSignal->LongSignalTimeComplex()[0][0], phaseLagRF, fPhaseLagLO, timeStamp);
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

    double fRF_frequency;
    double fLO_frequency;
    double fAcquisitionRate;
    double fFilterRate;
    double fAmplitude;
    double fPhaseLagLO;



};

bool parsePlaneWaveFIR(test_app& the_main)
{
	TestParameterHandler* p1 = TestParameterHandler::getInstance();
    CLI11_PARSE( the_main, p1->GetArgc(), p1->GetArgv() );
	return true;
}



TEST_CASE( "LMCPlaneWaveFIR with default parameter values (pass)", "[single-file]" )
{
	test_app the_main;
	if (!parsePlaneWaveFIR(the_main)) exit(-1);

	testLMCPlaneWaveFIR aTestLMCPlaneWaveFIR;
	if ( !aTestLMCPlaneWaveFIR.Configure() )
	{
		LWARN(testlog,"testLMCPlaneWaveFIR was not configured correctly.");
	    REQUIRE( 0 > 1 );
	}

    /* initialize time series */
    Signal* aSignal = new Signal();
    int N0 = aTestLMCPlaneWaveFIR.fTFReceiverHandler->GetFilterSize();
    aTestLMCPlaneWaveFIR.fFilterRate = (1./aTestLMCPlaneWaveFIR.fTFReceiverHandler->GetFilterResolution());
    aSignal->Initialize( N0 , 1 );

    double firGainMax = 0.;
    for (unsigned iHarmonic=1; iHarmonic<3; iHarmonic++)
    {
        for (unsigned rfStep=0; rfStep<25; rfStep++) // frequency sweep
        {
            aTestLMCPlaneWaveFIR.fRF_frequency = iHarmonic*20.9e9 + 0.5e9*rfStep;

	        double convolutionMag = 0.;
            for (unsigned i=0; i<1000; i++)  // time stamps
            {
               /* populate time series and convolve it with the FIR filter */
        	    double timeStamp = i/aTestLMCPlaneWaveFIR.fAcquisitionRate;
                aTestLMCPlaneWaveFIR.PopulateSignal(aSignal, N0, timeStamp);
                double convolution = aTestLMCPlaneWaveFIR.fTFReceiverHandler->ConvolveWithFIRFilter(aTestLMCPlaneWaveFIR.SignalToDeque(aSignal));
                if (fabs(convolution) > convolutionMag)
                {
        	        convolutionMag = convolution;
                }
            } // i
            // https://www.antenna-theory.com/definitions/antennafactor.php
            double firGain = 10.*log10(pow(1./(aTestLMCPlaneWaveFIR.fAmplitude/convolutionMag/9.73*(3.e8/aTestLMCPlaneWaveFIR.fRF_frequency)),2.));
            if (firGain > firGainMax) firGainMax = firGain;
			LPROG( testlog, "FIR gain at frequency " << aTestLMCPlaneWaveFIR.fRF_frequency << " is " << firGain << " dB" );
//            printf("firGain at frequency %g is %g dB\n", aTestLMCPlaneWaveFIR.fRF_frequency, firGain);
        } // rfStep
    }

    delete aSignal;

    REQUIRE( firGainMax > 0. );
}



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

class testCavity_app : public main_app
{
    public:
        testCavity_app() :
            main_app(),
			fDHOTimeResolution(1.e-8),
			fDHOThresholdFactor(0.01),
			fCavityFrequency(1.067e9),
			fCavityQ(1000.)
        {
            add_option("-r,--dho-time-resolution", fDHOTimeResolution, "Time resolution used in Green's function (s).");
            add_option("-m,--dho-threshold-factor", fDHOThresholdFactor, "Minimum fractional threshold of Green's function used to calculate FIR.");
            add_option("-f,--dho-cavity-frequency", fCavityFrequency, "Cavity resonant frequency (Hz).");
            add_option("-q,--dho-cavity-Q", fCavityQ, "Cavity Q.");
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


    private:
        double fDHOTimeResolution;
        double fDHOThresholdFactor;
        double fCavityFrequency;
        double fCavityQ;
};


class testLMCCavity
{
public:

	testLMCCavity():
        fRF_frequency(0.),
        fAcquisitionRate(0.),
		fFilterRate(0.),
		fTFReceiverHandler( 0 ),
		fAnalyticResponseFunction( 0 ),
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
		testCavity_app the_main;
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

        fRF_frequency = 0.; // Hz
        fAcquisitionRate = 201.e6; // Hz
		return true;
	}



	bool PopulateSignal(Signal* aSignal, int N0, double timeStamp)
	{

        double voltage_phase = 0.;
        double phaseLagRF = 2.*LMCConst::Pi() * timeStamp/(1./fRF_frequency);

        for( unsigned index = 0; index < N0; ++index )
        {
            voltage_phase = phaseLagRF + 2.*LMCConst::Pi()*fRF_frequency*(double)index/fFilterRate;

            aSignal->LongSignalTimeComplex()[index][0] = cos(voltage_phase);
            aSignal->LongSignalTimeComplex()[index][1] = cos(-LMCConst::Pi()/2. + voltage_phase);
        }
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
    scarab::param_node* param_0;


    double fRF_frequency;
    double fAcquisitionRate;
    double fFilterRate;

};


int parseCavity(testCavity_app& the_main)
{
	TestParameterHandler* p1 = TestParameterHandler::getInstance();
    CLI11_PARSE( the_main, p1->GetArgc(), p1->GetArgv() );
	return 0;
}


TEST_CASE( "testLMCCavity with default parameter values (pass)", "[single-file]" )
{
	testCavity_app the_main;
	parseCavity(the_main);

	testLMCCavity aTestLMCCavity;
	aTestLMCCavity.AddParam( "dho-time-resolution", the_main.GetDHOTimeResolution() );
	aTestLMCCavity.AddParam( "dho-threshold-factor", the_main.GetDHOThresholdFactor() );
	aTestLMCCavity.AddParam( "dho-cavity-frequency", the_main.GetCavityFrequency() );
	aTestLMCCavity.AddParam( "dho-cavity-Q", the_main.GetCavityQ() );

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

    double qInferred = 0.;
    double maxGain = 0.;
    double rfSpanSweep = 3. * the_main.GetCavityFrequency() / the_main.GetCavityQ();
    double rfStepSize = 0.00005 * the_main.GetCavityFrequency();
    int nSteps = rfSpanSweep / rfStepSize;

        for (int rfStep=-nSteps/2; rfStep<nSteps/2; rfStep++) // frequency sweep
        {
            aTestLMCCavity.fRF_frequency = the_main.GetCavityFrequency() + rfStepSize * rfStep;

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

            if (convolutionMag*convolutionMag > maxGain)
            {
            	maxGain = convolutionMag*convolutionMag;
            }
            else if ((convolutionMag*convolutionMag < 0.5*maxGain) && (qInferred == 0.))
            {
            	qInferred = the_main.GetCavityFrequency() /  (2.* rfStepSize * (rfStep-1));
            }
			LPROG( testlog, "Cavity GF gain at frequency " << aTestLMCCavity.fRF_frequency << " is " << convolutionMag );
        } // rfStep

    delete aSignal;

    LPROG( testlog, "Estimated Q is " << qInferred );
    LPROG( testlog, "Expected Q is " << the_main.GetCavityQ() );
    REQUIRE( qInferred > 0. );
}



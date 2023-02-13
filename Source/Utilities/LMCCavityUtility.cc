/*
 * LMCCavityUtility.cc
 *
 *  Created on: Feb 10, 2023
 *      Author: pslocum
 */

#include "LMCCavityUtility.hh"
#include "logger.hh"

using namespace scarab;

namespace locust
{
    LOGGER( testlog, "CavityUtility" );

    CavityUtility::CavityUtility():
    		fRF_frequency(0.),
			fFilterRate(0.),
			fTFReceiverHandler( 0 ),
			fAnalyticResponseFunction( 0 ),
			fparam_0( new param_node() )
    {
    }

    CavityUtility::~CavityUtility()
    {
    }

	bool CavityUtility::Configure()
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

		return true;
	}

	void CavityUtility::AddParam(std::string aString, double aValue)
	{
		if ( aValue != 0. )
		{
            fparam_0->add(aString.c_str(), aValue);
		}
	}

	const scarab::param_node* CavityUtility::GetParams()
	{
	    return fparam_0;
	}


	bool CavityUtility::PopulateSignal(Signal* aSignal, int N0)
	{

        double voltage_phase = 0.;

        for( unsigned index = 0; index < N0; ++index )
        {
            voltage_phase = 2.*LMCConst::Pi()*fRF_frequency*(double)index/fFilterRate;

            aSignal->LongSignalTimeComplex()[index][0] = cos(voltage_phase);
            aSignal->LongSignalTimeComplex()[index][1] = cos(-LMCConst::Pi()/2. + voltage_phase);
        }
        return true;
	}



	std::deque<double> CavityUtility::SignalToDeque(Signal* aSignal)
	{
	    std::deque<double> incidentSignal;
	    for (unsigned i=0; i<fTFReceiverHandler->GetFilterSize(); i++)
	    {
	    	incidentSignal.push_back(aSignal->LongSignalTimeComplex()[i][0]);
	    }
	    return incidentSignal;
	}


    bool CavityUtility::CheckCavityQ( double dhoTimeResolution, double dhoThresholdFactor, double dhoCavityFrequency, double dhoCavityQ)
    {
    	AddParam( "dho-time-resolution", dhoTimeResolution );
    	AddParam( "dho-threshold-factor", dhoThresholdFactor );
    	AddParam( "dho-cavity-frequency", dhoCavityFrequency );
    	AddParam( "dho-cavity-Q", dhoCavityQ );
    	if (!Configure())
    	{
    		LERROR(testlog,"Cavity was not configured correctly.");
    	    exit(-1);
    	}

        /* initialize time series */
        Signal* aSignal = new Signal();
        int N0 = fTFReceiverHandler->GetFilterSize();
        fFilterRate = (1./fTFReceiverHandler->GetFilterResolution());
        aSignal->Initialize( N0 , 1 );

        double qInferred = 0.;
        double maxGain = 0.;
        double rfSpanSweep = 3. * dhoCavityFrequency / dhoCavityQ;
        double rfStepSize = 0.00005 * dhoCavityFrequency;
        int nSteps = rfSpanSweep / rfStepSize;

            for (int rfStep=-nSteps/2; rfStep<nSteps/2; rfStep++) // frequency sweep
            {
                fRF_frequency = dhoCavityFrequency + rfStepSize * rfStep;

    	        double convolutionMag = 0.;
                for (unsigned i=0; i<1000; i++)  // time stamps
                {
                    // populate time series and convolve it with the FIR filter
                    PopulateSignal(aSignal, N0);
                	std::pair<double,double> convolutionPair = fTFReceiverHandler->ConvolveWithComplexFIRFilter(SignalToDeque(aSignal));
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
                	qInferred = dhoCavityFrequency /  (2.* rfStepSize * (rfStep-1));
                }
    			LDEBUG( testlog, "Cavity GF gain at frequency " << fRF_frequency << " is " << convolutionMag );
            } // rfStep

        delete aSignal;

        LPROG( testlog, "\nSummary:");
    	LPROG( testlog, "dho-threshold-factor is " << dhoThresholdFactor );
    	LPROG( testlog, "dho-time-resolution is " << dhoTimeResolution );
    	LPROG( testlog, "dho-cavity-frequency is " << dhoCavityFrequency );
    	LPROG( testlog, "dho-cavity-Q is " << dhoCavityQ );


        LPROG( testlog, "Estimated Q is " << qInferred );
        LPROG( testlog, "Expected Q is " << dhoCavityQ );

        if ( fabs( 1. - qInferred / dhoCavityQ ) < 0.05 )
        {
            LPROG( testlog, "The cavity Q has been configured correctly." );
        	return true;
        }
        else
        {
        	LERROR( testlog, "The Q value is " << qInferred << " but was supposed to be " << dhoCavityQ);
        	return false;
        }
    }


} /* namespace locust */

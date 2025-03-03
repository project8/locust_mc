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
			fExpandFactor( 1.0 ),
			fWriteOutputFile( false ),
			fTFReceiverHandler( 0 ),
			fAnalyticResponseFunction( 0 ),
			fparam_0( new param_node() ),
			fOutputPath( TOSTRING(PB_OUTPUT_DIR) )
    {
    }

    CavityUtility::~CavityUtility()
    {
    }

	bool CavityUtility::Configure(int bTE, int l, int m, int n)
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

		if ( !fTFReceiverHandler->ConvertAnalyticGFtoFIR({{bTE,l,m,n}}, fAnalyticResponseFunction->GetGFarray({{bTE,l,m,n}})) )
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

	void CavityUtility::SetExpandFactor(double aFactor)
	{
		fExpandFactor = aFactor;
	}

	void CavityUtility::SetOutputFile(bool aFlag)
	{
		fWriteOutputFile = aFlag;
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



	std::deque<double> CavityUtility::SignalToDeque(int bTE, int l, int m, int n, Signal* aSignal)
	{
	    std::deque<double> incidentSignal;
	    for (unsigned i=0; i<fTFReceiverHandler->GetFilterSizeArray(bTE,l,m,n); i++)
	    {
	    	incidentSignal.push_back(aSignal->LongSignalTimeComplex()[i][0]);
	    }
	    return incidentSignal;
	}


	void CavityUtility::SetOutputPath( std::string aPath )
	{
		fOutputPath = aPath;
	}

	bool CavityUtility::WriteRootHisto(int npoints, double* freqArray, double* gainArray)
	{
	#ifdef ROOT_FOUND
		char cBufferFileName[60];
		int n = sprintf(cBufferFileName, "%s/UnitTestOutput.root", fOutputPath.c_str());
		const char *cFileName = cBufferFileName;
		FileWriter* aRootHistoWriter = RootHistoWriter::get_instance();
		aRootHistoWriter->SetFilename(cFileName);
		aRootHistoWriter->OpenFile("RECREATE");
		TH1D* aHisto = new TH1D("cavityHisto", "Green's function; frequency (Hz); 10 log10(|A|^{2})", npoints, freqArray[0], freqArray[npoints-1]);
		for (unsigned i=0; i<npoints; i++)
		{
			aHisto->SetBinContent(i+1, 10.*log10(gainArray[i]));
		}
		aRootHistoWriter->Write1DHisto(aHisto);
		aRootHistoWriter->CloseFile();
	#endif
		return true;
	}


    bool CavityUtility::CheckCavityQ( int nModes, int bTE, int l, int m, int n, double dhoTimeResolution, double dhoThresholdFactor, double dhoCavityFrequency, double dhoCavityQ)
    {
    	AddParam( "dho-time-resolution", dhoTimeResolution );
    	AddParam( "dho-threshold-factor", dhoThresholdFactor );
    	AddParam( "dho-cavity-frequency", dhoCavityFrequency );
    	AddParam( "dho-cavity-Q", dhoCavityQ );
    	AddParam( "n-modes", nModes );
    	if (!Configure(bTE, l, m, n))
    	{
    		LERROR(testlog,"Cavity was not configured correctly.");
    	    exit(-1);
    	}

        /* initialize time series */
        Signal* aSignal = new Signal();
        int N0 = fTFReceiverHandler->GetFilterSizeArray(bTE, l, m, n);
        fFilterRate = (1./fTFReceiverHandler->GetFilterResolutionArray(bTE, l, m, n));
        aSignal->Initialize( N0 , 1 );

        double qInferred = 0.;
        double maxGain = 0.;
        double rfSpanSweep = 3. * dhoCavityFrequency / dhoCavityQ;
        double rfStepSize = 0.00005 * dhoCavityFrequency;
        int nSteps = fExpandFactor * rfSpanSweep / rfStepSize;
        double* freqArray = new double[nSteps];
        double* gainArray = new double[nSteps];

        for (int i=0; i<nSteps; i++) // frequency sweep
        {
        	int rfStep = -nSteps/2/fExpandFactor + i;
        	fRF_frequency = dhoCavityFrequency + rfStepSize * rfStep;
        	double convolutionMag = 0.;
        	// populate time series and convolve it with the FIR filter
        	PopulateSignal(aSignal, N0);
        	std::pair<double,double> convolutionPair = fTFReceiverHandler->ConvolveWithComplexFIRFilterArray(bTE, l, m, n,SignalToDeque(bTE, l, m, n, aSignal));

        	if (fabs(convolutionPair.first) > convolutionMag)
        	{
        	    convolutionMag = convolutionPair.first;
        	}

        	freqArray[i] = fRF_frequency;
        	gainArray[i] = convolutionMag*convolutionMag;

        	if (convolutionMag*convolutionMag > maxGain)
        	{
        		maxGain = convolutionMag*convolutionMag;
        		qInferred = 0.;
        	}
        	else if ((convolutionMag*convolutionMag < 0.5*maxGain) && (qInferred == 0.))
        	{
        		qInferred = dhoCavityFrequency /  (2.* rfStepSize * (rfStep-1));
        	}
        	LPROG( testlog, "Cavity GF gain at frequency " << fRF_frequency << " is " << convolutionMag );
        }

#ifdef ROOT_FOUND
        if (fWriteOutputFile) WriteRootHisto(nSteps, freqArray, gainArray);
#endif
        aSignal->Reset();
        delete[] freqArray;
        delete[] gainArray;
        delete fTFReceiverHandler;
        delete fAnalyticResponseFunction;

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

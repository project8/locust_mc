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

		LPROG(testlog, "Configuring DHO in cavity utitlity.");
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

	void CavityUtility::AddParam(std::string aString, std::string aValue)
	{
	    if (!aValue.empty())
	    {
	        fparam_0->add(aString.c_str(), aValue.c_str());
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



	std::deque<double> CavityUtility::SignalToDeque(int bTE, int l, int m, int n, Signal* aSignal, int startIdx, int endIdx)
	{
	    std::deque<double> incidentSignal;
	    for (unsigned i=startIdx; i<endIdx; i++)
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
        	std::pair<double,double> convolutionPair = fTFReceiverHandler->ConvolveWithComplexFIRFilterArray(bTE, l, m, n,SignalToDeque(bTE, l, m, n, aSignal, 0, N0));

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

	bool CavityUtility::CheckCavityQNorm( int nModes, int bTE, int l, int m, int n, std::string dhoTimeResolution, double dhoThresholdFactor, double dhoCavityFrequency, double dhoCavityQ)
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

		// Magic to convert power-norm-specific dho string to float
		double ddhoTimeResolution;
		std::string fraction = dhoTimeResolution;
		size_t delimiterPos = fraction.find('/');
		if (delimiterPos != std::string::npos) {
			double numerator = std::stod(fraction.substr(0, delimiterPos));
			double denominator = std::stod(fraction.substr(delimiterPos + 1));
			ddhoTimeResolution = numerator / denominator;
		}

        /* initialize time series */
        Signal* aSignal = new Signal();
        int N0 = 3 * dhoCavityQ/dhoCavityFrequency / ddhoTimeResolution; // Calculate E field for x ring-up times
        fFilterRate = (1./fTFReceiverHandler->GetFilterResolutionArray(bTE, l, m, n));
        aSignal->Initialize( N0 , 1 );

		LPROG( testlog, "(norm) N0 is " << N0 );

        double qInferred = 0.;
        double maxGain = 0.;
        double rfSpanSweep = 3. * dhoCavityFrequency / dhoCavityQ;
        double rfStepSize = 0.000005 * dhoCavityFrequency;
        int nSteps = fExpandFactor * rfSpanSweep / rfStepSize;
        double* freqArray = new double[nSteps];
        double* gainArray = new double[nSteps];

		// Size of GF in DHO is finite. May need to propagate in large steps
		int numGFBins = fTFReceiverHandler->GetFilterSizeArray(bTE, l, m, n) - 1; // -1 subtracts wprime, fBfactor bins
		int nLargeSteps = 1 + std::floor( (N0-1) / numGFBins);

        for (int i=0; i<nSteps; i++) // frequency sweep
        {
        	int rfStep = -nSteps/2/fExpandFactor + i;
        	fRF_frequency = dhoCavityFrequency + rfStepSize * rfStep;
        	double convolutionMag = 0.;
        	// populate time series and convolve it with the FIR filter
        	PopulateSignal(aSignal, N0);
			// Size of GF in DHO is finite. May need to propagate in large steps
			for (int j=0; j<nLargeSteps; j++)
			{
				int startIdx = j*numGFBins;
				int endIdx;
				if (std::floor((N0 - startIdx) / numGFBins)) // if there is another large step after this
				{
					endIdx = (j+1)*numGFBins;
				}
				else // if this is the last large step
				{
					endIdx = startIdx + (N0 % numGFBins);
				}
				double tProp = startIdx * fTFReceiverHandler->GetFilterResolutionArray(bTE, l, m, n);

				fTFReceiverHandler->ComputeFields(bTE, l, m, n,SignalToDeque(bTE, l, m, n, aSignal, startIdx, endIdx), tProp);
			}
			// Grab the last computed complex E field
        	const auto& convolutionPair = (fTFReceiverHandler->GetEfield()[bTE][l][m][n]).back();
			double eFieldAbs = sqrt( convolutionPair[0]*convolutionPair[0] + convolutionPair[1]*convolutionPair[1] );

        	if (eFieldAbs > convolutionMag)
        	{
        	    convolutionMag = eFieldAbs;
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
        	LPROG( testlog, "(norm) Cavity GF gain at frequency " << fRF_frequency << " is " << convolutionMag );

			// Reset the field in the cavity to zero.
			fTFReceiverHandler->SetLastElementToZero(bTE, l, m, n);
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
        LPROG( testlog, "(norm) dho-threshold-factor is " << dhoThresholdFactor );
        LPROG( testlog, "(norm) dho-time-resolution is " << ddhoTimeResolution );
        LPROG( testlog, "(norm) dho-cavity-frequency is " << dhoCavityFrequency );
        LPROG( testlog, "(norm) dho-cavity-Q is " << dhoCavityQ );
        LPROG( testlog, "(norm) Estimated Q is " << qInferred );
        LPROG( testlog, "(norm) Expected Q is " << dhoCavityQ );

        if ( fabs( 1. - qInferred / dhoCavityQ ) < 0.05 )
        {
            LPROG( testlog, "(norm) The cavity Q has been configured correctly." );
        	return true;
        }
        else
        {
        	LERROR( testlog, "(norm) The Q value is " << qInferred << " but was supposed to be " << dhoCavityQ);
        	return false;
        }
    }


} /* namespace locust */

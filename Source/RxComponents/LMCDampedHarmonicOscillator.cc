/*
 * LMCDampedHarmonicOscillator.cc
 *
 *  Created on: Jul 6, 2022
 *      Author: pslocum
 */

#include "LMCDampedHarmonicOscillator.hh"
#include "logger.hh"

namespace locust
{

	LOGGER( lmclog, "DampedHarmonicOscillator" );

    DampedHarmonicOscillator::DampedHarmonicOscillator():
    		fMaxNBins( 20000 ),
			fTimeResolution( {{{1.e-10}}} ),
			fNModes( 2 ),
			fCavityFrequency( {{{1.067e9}}} ),
			fCavityQ( {{{1000}}} ),
			fThresholdFactor ( {{{0.25}}} ),
			fCavityDampingFactor( {{{0.}}} ),
			fBFactor( {{{0.}}} ),
			fHannekePowerFactor( {{{1.}}} )
    {

	fCavityFrequency.resize(fNModes);
	fCavityQ.resize(fNModes);
	//fTimeResolution.resize(fNModes);
	fThresholdFactor.resize(fNModes);
	fBFactor.resize(fNModes);
	fCavityDampingFactor.resize(fNModes);
	fHannekePowerFactor.resize(fNModes);
	for(int l=0; l<fNModes; l++)
	{
		fCavityFrequency[l].resize(fNModes);
         	fCavityQ[l].resize(fNModes);
        	//fTimeResolution[l].resize(fNModes);
        	fThresholdFactor[l].resize(fNModes);
          	fBFactor[l].resize(fNModes);
		fCavityDampingFactor[l].resize(fNModes);
		fHannekePowerFactor[l].resize(fNModes);
		for(int m=0; m<fNModes; m++)
		{
			fCavityFrequency[l][m].resize(fNModes);
			fCavityQ[l][m].resize(fNModes);
        		//fTimeResolution[l][m].resize(fNModes);
        		fThresholdFactor[l][m].resize(fNModes);
			fBFactor[l][m].resize(fNModes);
			fCavityDampingFactor[l][m].resize(fNModes);
			fHannekePowerFactor[l][m].resize(fNModes);
			for(int n=0; n<fNModes; n++)
                  	{
				fCavityFrequency[l][m][n] = 1.067e9;
				fCavityQ[l][m][n] = 1000.;
        			//fTimeResolution[l][m][n] = 1.e-10;
        			fThresholdFactor[l][m][n] = 0.25;
				fBFactor[l][m][n] = 0.;
				fCavityDampingFactor[l][m][n] = 0.;
				fHannekePowerFactor[l][m][n] = 1.;
			}
		}	
	}

    }
    DampedHarmonicOscillator::~DampedHarmonicOscillator() {}

    bool DampedHarmonicOscillator::Configure( const scarab::param_node& aParam )
    {
    	if( !AnalyticResponseFunction::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring AnalyticResponseFunction class from DampedHarmonicOscillator subclass");
    	}
    	if( aParam.has( "n-modes" ) )
        {
        	fNModes = aParam["n-modes"]().as_int();
        }
    	if( aParam.has( "dho-max-nbins" ) )
    	{
    		fMaxNBins = aParam["dho-max-nbins"]().as_int();
    	}
    	if( aParam.has( "dho-time-resolution" ) )
        {
                fTimeResolution.resize(fNModes);
                for (unsigned l=0; l<fNModes; l++)
                {
                        fTimeResolution[l].resize(fNModes);
                        for (unsigned m=0; m<fNModes; m++)
                        {
                                fTimeResolution[l][m].resize(fNModes);
                                for (unsigned n=0; n<fNModes; n++)
                                {
                                        SetDHOTimeResolution( l, m, n, aParam["dho-time-resolution"]().as_double() );
                                }
                        }
        	}
	}     
    	if( aParam.has( "dho-threshold-factor" ) )
        {  
		fThresholdFactor.resize(fNModes); 
                for (unsigned l=0; l<fNModes; l++)
                {   
                        fThresholdFactor[l].resize(fNModes);
                        for (unsigned m=0; m<fNModes; m++)
                        {   
                                fThresholdFactor[l][m].resize(fNModes);
                                for (unsigned n=0; n<fNModes; n++)
                                {
					SetDHOThresholdFactor( l, m, n, aParam["dho-threshold-factor"]().as_double() );   
                                }   
                        }   
                }   
        }   
    	if( aParam.has( "dho-cavity-frequency" ) )
    	{
		fCavityFrequency.resize(fNModes);
		for (unsigned l=0; l<fNModes; l++)
		{
			fCavityFrequency[l].resize(fNModes);
			for (unsigned m=0; m<fNModes; m++)
			{
				fCavityFrequency[l][m].resize(fNModes);
				for (unsigned n=0; n<fNModes; n++)
				{
					SetCavityFrequency( l, m, n, aParam["dho-cavity-frequency"]().as_double() );
				}
			}
		}
    	}
	if( aParam.has( "dho-cavity-frequency-TE011" ) ) 
        {
		SetCavityFrequency( 0, 1, 1, aParam["dho-cavity-frequency-TE011"]().as_double() );	
	}
        if( aParam.has( "dho-cavity-frequency-TE111" ) )
        {       
                SetCavityFrequency( 1, 1, 1, aParam["dho-cavity-frequency-TE111"]().as_double() );
        }
    	if( aParam.has( "dho-cavity-Q" ) )
        {   
                fCavityQ.resize(fNModes);
                for (unsigned l=0; l<fNModes; l++)
                {   
                        fCavityQ[l].resize(fNModes);
                        for (unsigned m=0; m<fNModes; m++)
                        {   
                                fCavityQ[l][m].resize(fNModes);
                                for (unsigned n=0; n<fNModes; n++)
                                {   
                                        SetCavityQ( l, m, n, aParam["dho-cavity-Q"]().as_double() );
                                }   
                        }   
                }  
        }   
    	fTFReceiverHandler = new TFReceiverHandler;
    	if(!fTFReceiverHandler->Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver FIRHandler class");
    		exit(-1);
    		return false;
    	}
    	if ( !Initialize() ) return false;
    	else return true;
    }

    bool DampedHarmonicOscillator::Initialize()
    {
	fCavityOmega.resize(fNModes);
	fCavityDampingFactor.resize(fNModes);
	fBFactor.resize(fNModes);
        fCavityOmegaPrime.resize(fNModes);
        for (unsigned l=0; l<fNModes; l++)
        {   
        	fCavityOmega[l].resize(fNModes);
        	fCavityDampingFactor[l].resize(fNModes);
        	fBFactor[l].resize(fNModes);
        	fCavityOmegaPrime[l].resize(fNModes);
                for (unsigned m=0; m<fNModes; m++)
                {   
        		fCavityOmega[l][m].resize(fNModes);
        		fCavityDampingFactor[l][m].resize(fNModes);
        		fBFactor[l][m].resize(fNModes);
        		fCavityOmegaPrime[l][m].resize(fNModes); 
                        for (unsigned n=0; n<fNModes; n++)
                        {   
        			fCavityOmega[l][m][n] = fCavityFrequency[l][m][n] * 2. * LMCConst::Pi();
        			fCavityDampingFactor[l][m][n] = 1. / 2. / fCavityQ[l][m][n];
        			fBFactor[l][m][n] = fCavityDampingFactor[l][m][n] * fCavityOmega[l][m][n];
        			fCavityOmegaPrime[l][m][n] = sqrt( fCavityOmega[l][m][n]*fCavityOmega[l][m][n] - fBFactor[l][m][n]*fBFactor[l][m][n] );
				if (!GenerateGreensFunction(l,m,n)) return false;
                        }   
                }   
        }   
    	return true;
    }

    void DampedHarmonicOscillator::SetCavityQ( int l, int m, int n, double aQ )
    {
    	fCavityQ[l][m][n] = aQ;
    }
    double DampedHarmonicOscillator::GetCavityQ(int l, int m, int n)
    {
    	return fCavityQ[l][m][n];
    }
    void DampedHarmonicOscillator::SetCavityFrequency( int l, int m, int n, double aFrequency )
    {
    	fCavityFrequency[l][m][n] = aFrequency;
    }
    double DampedHarmonicOscillator::GetCavityFrequency(int l, int m, int n)
    {
    	return fCavityFrequency[l][m][n];
    }
    void DampedHarmonicOscillator::SetDHOTimeResolution(int l, int m, int n, double aTimeResolution )
    {
    	fTimeResolution[l][m][n] = aTimeResolution;
    }
    double DampedHarmonicOscillator::GetDHOTimeResolution(int l, int m, int n)
    {
    	return fTimeResolution[l][m][n];
    }
    void DampedHarmonicOscillator::SetDHOThresholdFactor(int l, int m, int n, double aThresholdFactor )
    {
    	fThresholdFactor[l][m][n] = aThresholdFactor;
    }
    double DampedHarmonicOscillator::GetDHOThresholdFactor(int l, int m, int n)
    {
    	return fThresholdFactor[l][m][n];
    }



    double DampedHarmonicOscillator::ExpDecayTerm(int l, int m, int n, double t)
    {
    	double ExpDecayTerm = exp( -fBFactor[l][m][n] * t);
    	return ExpDecayTerm;
    }

    std::pair<double,double> DampedHarmonicOscillator::GreensFunction(int l, int m, int n, double t)
    {
    	//double GreensFunctionValueReal = ExpDecayTerm(t) * sin( fCavityOmegaPrime * t) / fCavityOmegaPrime;
    	//double GreensFunctionValueImag = -ExpDecayTerm(t) * cos( fCavityOmegaPrime * t) / fCavityOmegaPrime;

    	// Modify Green's function for nominal gain of unity, keeping phase information unchanged.
    	// Power model could possibly be implemented as in here:

    	double GreensFunctionValueReal = fTimeResolution[l][m][n] * fHannekePowerFactor[l][m][n] * ExpDecayTerm(l, m, n, t) * sin( fCavityOmegaPrime[l][m][n] * t);
    	double GreensFunctionValueImag = -1. * fTimeResolution[l][m][n] * fHannekePowerFactor[l][m][n] * ExpDecayTerm(l, m, n, t) * cos( fCavityOmegaPrime[l][m][n] * t);
	
    	return std::make_pair(GreensFunctionValueReal,GreensFunctionValueImag);
    }


    double DampedHarmonicOscillator::NormFactor(int l, int m, int n, double aDriveFrequency)
    {
		if (!fTFReceiverHandler->ConvertAnalyticGFtoFIR(l, m, n, GetGFarray(l, m, n)))
		{
			LERROR(lmclog,"GF->FIR was not generated in DHO::NormFactor.");
			exit(-1);
			return false;
		}
        /* initialize time series */
        Signal* aSignal = new Signal();
        int N0 = GetGFarray(l, m, n).size();
        aSignal->Initialize( N0 , 1 );
        double convolutionMag = 0.;

        //for (unsigned i=0; i<1000; i++)  // time stamps
        for (unsigned i=0; i<1; i++) 
        {
            // populate time series and convolve it with the FIR filter
            PopulateCalibrationSignal(l, m, n, aSignal, N0, aDriveFrequency);
        	std::pair<double,double> convolutionPair = fTFReceiverHandler->ConvolveWithComplexFIRFilterArray(l, m, n, SignalToDeque(aSignal));
            if (fabs(convolutionPair.first) > convolutionMag)
            {
    	        convolutionMag = convolutionPair.first;
            }
        } //

        delete aSignal;
        return convolutionMag;

    }

	bool DampedHarmonicOscillator::PopulateCalibrationSignal(int l, int m, int n, Signal* aSignal, int N0, double aDriveFrequency)
	{

        double voltage_phase = 0.;

        for( unsigned index = 0; index < N0; ++index )
        {
            voltage_phase = 2.*LMCConst::Pi()*aDriveFrequency*(double)index*fTimeResolution[l][m][n];

            aSignal->LongSignalTimeComplex()[index][0] = cos(voltage_phase);
            aSignal->LongSignalTimeComplex()[index][1] = cos(-LMCConst::Pi()/2. + voltage_phase);
        }
        return true;
	}

	std::deque<double> DampedHarmonicOscillator::SignalToDeque(Signal* aSignal)
	{
	    std::deque<double> incidentSignal;
	    for (unsigned i=0; i<fTFReceiverHandler->GetFilterSize(); i++)
	    {
	    	incidentSignal.push_back(aSignal->LongSignalTimeComplex()[i][0]);
	    }
	    return incidentSignal;
	}


    bool DampedHarmonicOscillator::GenerateGreensFunction(int l, int m, int n)
    {
        std::vector<std::pair<double,std::pair<double,double> > > tGFArray;

    	int sizeCounter = 0;

    	for (unsigned i=0; i<fMaxNBins; i++)
    	{
    		double tValue = i * fTimeResolution[l][m][n];
    		tGFArray.push_back(std::make_pair(fTimeResolution[l][m][n],GreensFunction(l, m, n, tValue)));
    		sizeCounter += 1;
    		if ( ExpDecayTerm(l, m, n, tValue) < fThresholdFactor[l][m][n] * ExpDecayTerm(l, m, n, 0.) )
    		{
    			break;
    		}
    	}
	
    	//tGFArray.resize( sizeCounter );
    	std::reverse( tGFArray.begin(), tGFArray.end() );
    	SetGFarray(l, m, n, tGFArray ); // unnormalized.
    	double aNormFactor = NormFactor(l, m, n, fCavityFrequency[l][m][n]);
    	for (unsigned i=0; i<sizeCounter; i++)
    	{
    		tGFArray[i].second.first /= aNormFactor;
    		tGFArray[i].second.second /= aNormFactor;
    	}
    	SetGFarray(l, m, n, tGFArray); // now normalized.
    	if ( tGFArray.size() < 1 ) return false;
    	else return true;

    }




} /* namespace locust */


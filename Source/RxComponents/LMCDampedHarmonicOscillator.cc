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
			fTimeResolution( {{{{1.e-10}}}} ),
			fNModes( 2 ),
			fCavityFrequency( {{{{1.067e9}}}} ),
			fCavityQ( {{{{1000}}}} ),
			fThresholdFactor ( {{{{0.25}}}} ),
			fCavityDampingFactor( {{{{0.}}}} ),
			fBFactor( {{{{0.}}}} ),
			fHannekePowerFactor( {{{{1.}}}} )
    {

	
        fCavityFrequency.resize(2);
        fCavityQ.resize(2);
        fTimeResolution.resize(2);
        fThresholdFactor.resize(2);
        fBFactor.resize(2);
        fCavityDampingFactor.resize(2);
        fHannekePowerFactor.resize(2);
	for(int bTE=0; bTE<2; bTE++)
	{

		fCavityFrequency[bTE].resize(fNModes);
		fCavityQ[bTE].resize(fNModes);
		fTimeResolution[bTE].resize(fNModes);
		fThresholdFactor[bTE].resize(fNModes);
		fBFactor[bTE].resize(fNModes);
		fCavityDampingFactor[bTE].resize(fNModes);
		fHannekePowerFactor[bTE].resize(fNModes);
		for(int l=0; l<fNModes; l++)
		{
			fCavityFrequency[bTE][l].resize(fNModes);
         		fCavityQ[bTE][l].resize(fNModes);
        		fTimeResolution[bTE][l].resize(fNModes);
        		fThresholdFactor[bTE][l].resize(fNModes);
          		fBFactor[bTE][l].resize(fNModes);
			fCavityDampingFactor[bTE][l].resize(fNModes);
			fHannekePowerFactor[bTE][l].resize(fNModes);
			for(int m=0; m<fNModes; m++)
			{
				fCavityFrequency[bTE][l][m].resize(fNModes);
				fCavityQ[bTE][l][m].resize(fNModes);
        			fTimeResolution[bTE][l][m].resize(fNModes);
        			fThresholdFactor[bTE][l][m].resize(fNModes);
				fBFactor[bTE][l][m].resize(fNModes);
				fCavityDampingFactor[bTE][l][m].resize(fNModes);
				fHannekePowerFactor[bTE][l][m].resize(fNModes);
				for(int n=0; n<fNModes; n++)
                  		{
					fCavityFrequency[bTE][l][m][n] = 1.067e9;
					fCavityQ[bTE][l][m][n] = 1000.;
        				fTimeResolution[bTE][l][m][n] = 1.e-10;
        				fThresholdFactor[bTE][l][m][n] = 0.25;
					fBFactor[bTE][l][m][n] = 0.;
					fCavityDampingFactor[bTE][l][m][n] = 0.;
					fHannekePowerFactor[bTE][l][m][n] = 1.;
				}
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
		fTimeResolution.resize(2);
		for (unsigned bTE=0; bTE<2; bTE++)
		{
                	fTimeResolution[bTE].resize(fNModes);
                	for (unsigned l=0; l<fNModes; l++)
                	{
                        	fTimeResolution[bTE][l].resize(fNModes);
                        	for (unsigned m=0; m<fNModes; m++)
                        	{
                                	fTimeResolution[bTE][l][m].resize(fNModes);
                                	for (unsigned n=0; n<fNModes; n++)
                                	{
                                        	SetDHOTimeResolution(bTE, l, m, n, aParam["dho-time-resolution"]().as_double() );
					}
                                }
                        }
        	}
	}     
    	if( aParam.has( "dho-threshold-factor" ) )
        {  
		fThresholdFactor.resize(2);
		for (unsigned bTE=0; bTE<2; bTE++)
                {
			fThresholdFactor[bTE].resize(fNModes); 
                	for (unsigned l=0; l<fNModes; l++)
                	{   
                        	fThresholdFactor[bTE][l].resize(fNModes);
                        	for (unsigned m=0; m<fNModes; m++)
                        	{   
                                	fThresholdFactor[bTE][l][m].resize(fNModes);
                                	for (unsigned n=0; n<fNModes; n++)
                                	{
						SetDHOThresholdFactor(bTE, l, m, n, aParam["dho-threshold-factor"]().as_double() );   
					}
                                }   
                        }   
                }   
        }   
    	if( aParam.has( "dho-cavity-frequency" ) )
    	{
		fCavityFrequency.resize(2);
		for (unsigned bTE=0; bTE<2; bTE++)
                {
			fCavityFrequency[bTE].resize(fNModes);
			for (unsigned l=0; l<fNModes; l++)
			{
				fCavityFrequency[bTE][l].resize(fNModes);
				for (unsigned m=0; m<fNModes; m++)
				{
					fCavityFrequency[bTE][l][m].resize(fNModes);
					for (unsigned n=0; n<fNModes; n++)
					{
						SetCavityFrequency(bTE, l, m, n, aParam["dho-cavity-frequency"]().as_double() );
					}
				}
			}
		}
    	}
	if( aParam.has( "dho-cavity-frequency-TE011" ) ) 
        {
		SetCavityFrequency( 1, 0, 1, 1, aParam["dho-cavity-frequency-TE011"]().as_double() );	
	}
        if( aParam.has( "dho-cavity-frequency-TM111" ) )
        {       
                SetCavityFrequency( 0, 1, 1, 1, aParam["dho-cavity-frequency-TM111"]().as_double() );
        }
    	if( aParam.has( "dho-cavity-Q" ) )
        {   
		fCavityQ.resize(2);
		for (unsigned bTE=0; bTE<2; bTE++)
		{
                	fCavityQ[bTE].resize(fNModes);
                	for (unsigned l=0; l<fNModes; l++)
                	{   
                        	fCavityQ[bTE][l].resize(fNModes);
                        	for (unsigned m=0; m<fNModes; m++)
                        	{   
                                	fCavityQ[bTE][l][m].resize(fNModes);
                                	for (unsigned n=0; n<fNModes; n++)
                                	{   
                                        	SetCavityQ(bTE, l, m, n, aParam["dho-cavity-Q"]().as_double() );
					}
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
        fCavityOmega.resize(2);
        fCavityDampingFactor.resize(2);
        fBFactor.resize(2);
        fCavityOmegaPrime.resize(2);	

	for (unsigned bTE=0; bTE<2; bTE++)
	{
		fCavityOmega[bTE].resize(fNModes);
		fCavityDampingFactor[bTE].resize(fNModes);
		fBFactor[bTE].resize(fNModes);
        	fCavityOmegaPrime[bTE].resize(fNModes);
        	for (unsigned l=0; l<fNModes; l++)
        	{   
        		fCavityOmega[bTE][l].resize(fNModes);
        		fCavityDampingFactor[bTE][l].resize(fNModes);
        		fBFactor[bTE][l].resize(fNModes);
        		fCavityOmegaPrime[bTE][l].resize(fNModes);
                	for (unsigned m=0; m<fNModes; m++)
                	{   
        			fCavityOmega[bTE][l][m].resize(fNModes);
        			fCavityDampingFactor[bTE][l][m].resize(fNModes);
        			fBFactor[bTE][l][m].resize(fNModes);
        			fCavityOmegaPrime[bTE][l][m].resize(fNModes); 
                        	for (unsigned n=0; n<fNModes; n++)
                        	{   
        				fCavityOmega[bTE][l][m][n] = fCavityFrequency[bTE][l][m][n] * 2. * LMCConst::Pi();
        				fCavityDampingFactor[bTE][l][m][n] = 1. / 2. / fCavityQ[bTE][l][m][n];
        				fBFactor[bTE][l][m][n] = fCavityDampingFactor[bTE][l][m][n] * fCavityOmega[bTE][l][m][n];
        				fCavityOmegaPrime[bTE][l][m][n] = sqrt( fCavityOmega[bTE][l][m][n]*fCavityOmega[bTE][l][m][n] - fBFactor[bTE][l][m][n]*fBFactor[bTE][l][m][n] );
					if (!GenerateGreensFunction(bTE,l,m,n)) return false;
				}
                        }   
                }   
        }   
    	return true;
    }

    void DampedHarmonicOscillator::SetCavityQ(int bTE, int l, int m, int n, double aQ )
    {
    	fCavityQ[bTE][l][m][n] = aQ;
    }
    double DampedHarmonicOscillator::GetCavityQ(int bTE, int l, int m, int n)
    {
    	return fCavityQ[bTE][l][m][n];
    }
    void DampedHarmonicOscillator::SetCavityFrequency(int bTE, int l, int m, int n, double aFrequency )
    {
    	fCavityFrequency[bTE][l][m][n] = aFrequency;
    }
    double DampedHarmonicOscillator::GetCavityFrequency(int bTE, int l, int m, int n)
    {
    	return fCavityFrequency[bTE][l][m][n];
    }
    void DampedHarmonicOscillator::SetDHOTimeResolution(int bTE, int l, int m, int n, double aTimeResolution )
    {
    	fTimeResolution[bTE][l][m][n] = aTimeResolution;
    }
    double DampedHarmonicOscillator::GetDHOTimeResolution(int bTE, int l, int m, int n)
    {
    	return fTimeResolution[bTE][l][m][n];
    }
    void DampedHarmonicOscillator::SetDHOThresholdFactor(int bTE, int l, int m, int n, double aThresholdFactor )
    {
    	fThresholdFactor[bTE][l][m][n] = aThresholdFactor;
    }
    double DampedHarmonicOscillator::GetDHOThresholdFactor(int bTE, int l, int m, int n)
    {
    	return fThresholdFactor[bTE][l][m][n];
    }



    double DampedHarmonicOscillator::ExpDecayTerm(int bTE, int l, int m, int n, double t)
    {
    	double ExpDecayTerm = exp( -fBFactor[bTE][l][m][n] * t);
    	return ExpDecayTerm;
    }

    std::pair<double,double> DampedHarmonicOscillator::GreensFunction(int bTE, int l, int m, int n, double t)
    {
    	//double GreensFunctionValueReal = ExpDecayTerm(t) * sin( fCavityOmegaPrime * t) / fCavityOmegaPrime;
    	//double GreensFunctionValueImag = -ExpDecayTerm(t) * cos( fCavityOmegaPrime * t) / fCavityOmegaPrime;

    	// Modify Green's function for nominal gain of unity, keeping phase information unchanged.
    	// Power model could possibly be implemented as in here:

    	double GreensFunctionValueReal = fTimeResolution[bTE][l][m][n] * fHannekePowerFactor[bTE][l][m][n] * ExpDecayTerm(bTE, l, m, n, t) * sin( fCavityOmegaPrime[bTE][l][m][n] * t);
    	double GreensFunctionValueImag = -1. * fTimeResolution[bTE][l][m][n] * fHannekePowerFactor[bTE][l][m][n] * ExpDecayTerm(bTE, l, m, n, t) * cos( fCavityOmegaPrime[bTE][l][m][n] * t);
	
    	return std::make_pair(GreensFunctionValueReal,GreensFunctionValueImag);
    }


    double DampedHarmonicOscillator::NormFactor(int bTE, int l, int m, int n, double aDriveFrequency)
    {
		if (!fTFReceiverHandler->ConvertAnalyticGFtoFIR(bTE, l, m, n, GetGFarray(bTE, l, m, n)))
		{
			LERROR(lmclog,"GF->FIR was not generated in DHO::NormFactor.");
			exit(-1);
			return false;
		}
        /* initialize time series */
        Signal* aSignal = new Signal();
        int N0 = GetGFarray(bTE, l, m, n).size();
        aSignal->Initialize( N0 , 1 );
        double convolutionMag = 0.;

        //for (unsigned i=0; i<1000; i++)  // time stamps
        for (unsigned i=0; i<1; i++) 
        {
            // populate time series and convolve it with the FIR filter
            PopulateCalibrationSignal(bTE, l, m, n, aSignal, N0, aDriveFrequency);
        	std::pair<double,double> convolutionPair = fTFReceiverHandler->ConvolveWithComplexFIRFilterArray(bTE, l, m, n, SignalToDeque(aSignal));
            if (fabs(convolutionPair.first) > convolutionMag)
            {
    	        convolutionMag = convolutionPair.first;
            }
        } //

        delete aSignal;
        return convolutionMag;

    }

	bool DampedHarmonicOscillator::PopulateCalibrationSignal(int bTE, int l, int m, int n, Signal* aSignal, int N0, double aDriveFrequency)
	{

        double voltage_phase = 0.;

        for( unsigned index = 0; index < N0; ++index )
        {
            voltage_phase = 2.*LMCConst::Pi()*aDriveFrequency*(double)index*fTimeResolution[bTE][l][m][n];

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


    bool DampedHarmonicOscillator::GenerateGreensFunction(int bTE, int l, int m, int n)
    {
        std::vector<std::pair<double,std::pair<double,double> > > tGFArray;

    	int sizeCounter = 0;

    	for (unsigned i=0; i<fMaxNBins; i++)
    	{
    		double tValue = i * fTimeResolution[bTE][l][m][n];
    		tGFArray.push_back(std::make_pair(fTimeResolution[bTE][l][m][n],GreensFunction(bTE, l, m, n, tValue)));
    		sizeCounter += 1;
    		if ( ExpDecayTerm(bTE, l, m, n, tValue) < fThresholdFactor[bTE][l][m][n] * ExpDecayTerm(bTE, l, m, n, 0.) )
    		{
    			break;
    		}
    	}
	
    	//tGFArray.resize( sizeCounter );
    	std::reverse( tGFArray.begin(), tGFArray.end() );
    	SetGFarray(bTE, l, m, n, tGFArray ); // unnormalized.
    	double aNormFactor = NormFactor(bTE, l, m, n, fCavityFrequency[bTE][l][m][n]);
    	for (unsigned i=0; i<sizeCounter; i++)
    	{
    		tGFArray[i].second.first /= aNormFactor;
    		tGFArray[i].second.second /= aNormFactor;
    	}
    	SetGFarray(bTE, l, m, n, tGFArray); // now normalized.
    	if ( tGFArray.size() < 1 ) return false;
    	else return true;

    }


} /* namespace locust */


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
			fTimeResolution( 1.e-10 ),
			fCavityFrequency( 1.067e9 ),
			fCavityQ( 1000 ),
			fThresholdFactor ( 0.25 ),
			fCavityDampingFactor( 0. ),
			fBFactor( 0. ),
			fHannekePowerFactor( 1. )
    {}
    DampedHarmonicOscillator::~DampedHarmonicOscillator() {}

    bool DampedHarmonicOscillator::Configure( const scarab::param_node& aParam )
    {

    	if( !AnalyticResponseFunction::Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring AnalyticResponseFunction class from DampedHarmonicOscillator subclass");
    	}

    	if( aParam.has( "dho-max-nbins" ) )
    	{
    		fMaxNBins = aParam["dho-max-nbins"]().as_int();
    	}
    	if( aParam.has( "dho-time-resolution" ) )
    	{
    		SetDHOTimeResolution( aParam["dho-time-resolution"]().as_double() );
    	}
    	if( aParam.has( "dho-threshold-factor" ) )
    	{
    		SetDHOThresholdFactor( aParam["dho-threshold-factor"]().as_double() );
    	}
    	if( aParam.has( "dho-cavity-frequency" ) )
    	{
    		SetCavityFrequency( aParam["dho-cavity-frequency"]().as_double() );
    	}
    	if( aParam.has( "dho-cavity-Q" ) )
    	{
    		SetCavityQ( aParam["dho-cavity-Q"]().as_double() );
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
    	fCavityOmega = fCavityFrequency * 2. * LMCConst::Pi();
    	fCavityDampingFactor = 1. / 2. / fCavityQ;
    	fBFactor = fCavityDampingFactor * fCavityOmega;
    	fCavityOmegaPrime = sqrt( fCavityOmega*fCavityOmega - fBFactor*fBFactor );

    	if (!GenerateGreensFunction()) return false;
    	else return true;
    }

    void DampedHarmonicOscillator::SetCavityQ( double aQ )
    {
    	fCavityQ = aQ;
    }
    double DampedHarmonicOscillator::GetCavityQ()
    {
    	return fCavityQ;
    }
    void DampedHarmonicOscillator::SetCavityFrequency( double aFrequency )
    {
    	fCavityFrequency = aFrequency;
    }
    double DampedHarmonicOscillator::GetCavityFrequency()
    {
    	return fCavityFrequency;
    }
    void DampedHarmonicOscillator::SetDHOTimeResolution( double aTimeResolution )
    {
    	fTimeResolution = aTimeResolution;
    }
    double DampedHarmonicOscillator::GetDHOTimeResolution()
    {
    	return fTimeResolution;
    }
    void DampedHarmonicOscillator::SetDHOThresholdFactor( double aThresholdFactor )
    {
    	fThresholdFactor = aThresholdFactor;
    }
    double DampedHarmonicOscillator::GetDHOThresholdFactor()
    {
    	return fThresholdFactor;
    }



    double DampedHarmonicOscillator::ExpDecayTerm(double t)
    {
    	double ExpDecayTerm = exp( -fBFactor * t);
    	return ExpDecayTerm;
    }

    std::pair<double,double> DampedHarmonicOscillator::GreensFunction(double t)
    {
    	//double GreensFunctionValueReal = ExpDecayTerm(t) * sin( fCavityOmegaPrime * t) / fCavityOmegaPrime;
    	//double GreensFunctionValueImag = -ExpDecayTerm(t) * cos( fCavityOmegaPrime * t) / fCavityOmegaPrime;

    	// Modify Green's function for nominal gain of unity, keeping phase information unchanged.
    	// Power model could possibly be implemented as in here:

    	double GreensFunctionValueReal = fTimeResolution * fHannekePowerFactor * ExpDecayTerm(t) * sin( fCavityOmegaPrime * t);
    	double GreensFunctionValueImag = -1. * fTimeResolution * fHannekePowerFactor * ExpDecayTerm(t) * cos( fCavityOmegaPrime * t);

    	return std::make_pair(GreensFunctionValueReal,GreensFunctionValueImag);
    }


    double DampedHarmonicOscillator::NormFactor(double aDriveFrequency)
    {

		if (!fTFReceiverHandler->ConvertAnalyticGFtoFIR(1,0,1,1, GetGFarray()))
		{
			LERROR(lmclog,"GF->FIR was not generated in DHO::NormFactor.");
			exit(-1);
			return false;
		}

        /* initialize time series */
        Signal* aSignal = new Signal();
        int N0 = GetGFarray().size();
        aSignal->Initialize( N0 , 1 );
        double convolutionMag = 0.;

        for (unsigned i=0; i<1000; i++)  // time stamps
        {
            // populate time series and convolve it with the FIR filter
            PopulateCalibrationSignal(aSignal, N0, aDriveFrequency);

        	std::pair<double,double> convolutionPair = fTFReceiverHandler->ConvolveWithComplexFIRFilterArray(1,0,1,1,SignalToDeque(aSignal));

            if (fabs(convolutionPair.first) > convolutionMag)
            {
    	        convolutionMag = convolutionPair.first;
            }
        } // i

        delete aSignal;

        return convolutionMag;

    }

	bool DampedHarmonicOscillator::PopulateCalibrationSignal(Signal* aSignal, int N0, double aDriveFrequency)
	{

        double voltage_phase = 0.;

        for( unsigned index = 0; index < N0; ++index )
        {
            voltage_phase = 2.*LMCConst::Pi()*aDriveFrequency*(double)index*fTimeResolution;

            aSignal->LongSignalTimeComplex()[index][0] = cos(voltage_phase);
            aSignal->LongSignalTimeComplex()[index][1] = cos(-LMCConst::Pi()/2. + voltage_phase);
        }
        return true;
	}

	std::deque<double> DampedHarmonicOscillator::SignalToDeque(Signal* aSignal)
	{
	    std::deque<double> incidentSignal;
	    for (unsigned i=0; i<fTFReceiverHandler->GetFilterSizeArray(1,0,1,1); i++)
	    {
	    	incidentSignal.push_back(aSignal->LongSignalTimeComplex()[i][0]);
	    }
	    return incidentSignal;
	}


    bool DampedHarmonicOscillator::GenerateGreensFunction()
    {

        std::vector<std::pair<double,std::pair<double,double> > > tGFArray;

    	int sizeCounter = 0;

    	for (unsigned i=0; i<fMaxNBins; i++)
    	{
    		double tValue = i * fTimeResolution;
    		tGFArray.push_back(std::make_pair(fTimeResolution,GreensFunction(tValue)));
    		sizeCounter += 1;
    		if ( ExpDecayTerm(tValue) < fThresholdFactor * ExpDecayTerm(0.) )
    		{
    			break;
    		}
    	}

    	tGFArray.resize( sizeCounter );
    	std::reverse( tGFArray.begin(), tGFArray.end() );
    	SetGFarray( tGFArray ); // unnormalized.

    	double aNormFactor = NormFactor(fCavityFrequency);
    	for (unsigned i=0; i<sizeCounter; i++)
    	{
    		tGFArray[i].second.first /= aNormFactor;
    		tGFArray[i].second.second /= aNormFactor;
    	}

    	SetGFarray(tGFArray); // now normalized.

    	if ( tGFArray.size() < 1 ) return false;
    	else return true;

    }




} /* namespace locust */


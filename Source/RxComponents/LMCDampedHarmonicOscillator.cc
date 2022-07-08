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
    		fMaxNBins( 4000 ),
			fTimeResolution( 1.e-12 ),
			fCavityFrequency( 1.0e9 ),
			fCavityQ( 1000 ),
			fThresholdFactor ( 0.25 ),
			fCavityDampingFactor( 0. ),
			fBFactor( 0. )
    {}
    DampedHarmonicOscillator::~DampedHarmonicOscillator() {}

    bool DampedHarmonicOscillator::Configure( const scarab::param_node& aParam )
    {
    	if( aParam.has( "dho-max-nbins" ) )
    	{
    		fMaxNBins = aParam["dho-max-nbins"]().as_int();
    	}
    	if( aParam.has( "dho-time-resolution" ) )
    	{
    		fTimeResolution = aParam["dho-time-resolution"]().as_double();
    	}
    	if( aParam.has( "dho-cavity-frequency" ) )
    	{
    		fCavityFrequency = aParam["dho-cavity-frequency"]().as_double();
    	}
    	if( aParam.has( "dho-cavity-Q" ) )
    	{
    		fCavityQ = aParam["dho-cavity-Q"]().as_double();
    	}

    	Initialize();
    	return true;
    }

    bool DampedHarmonicOscillator::Initialize()
    {
    	fCavityOmega = fCavityFrequency * 2. * LMCConst::Pi();
    	fCavityDampingFactor = 1. / 2. / fCavityQ;
    	fBFactor = fCavityDampingFactor * fCavityOmega;
    	fCavityOmegaPrime = sqrt( fCavityOmega*fCavityOmega - fBFactor*fBFactor );

    	return true;
    }

    double DampedHarmonicOscillator::ExpDecayTerm(double t)
    {
    	double ExpDecayTerm = exp( -fBFactor * t);
    	return ExpDecayTerm;
    }

    double DampedHarmonicOscillator::GreensFunction(double t)
    {
    	double GreensFunctionValue = ExpDecayTerm(t) * sin( fCavityOmegaPrime * t) / fCavityOmegaPrime;
    	return GreensFunctionValue;
    }

    void DampedHarmonicOscillator::GenerateFIR()
    {

    	int sizeCounter = 0;
    	gfArray[0].first = fTimeResolution;
    	for (unsigned i=0; i<fMaxNBins; i++)
    	{
    		double tValue = i * fTimeResolution;
    		gfArray[i].second = GreensFunction( tValue );
    		sizeCounter += 1;
    		if ( ExpDecayTerm(tValue) < fThresholdFactor * ExpDecayTerm(0.) )
    		{
    			break;
    		}
    	}
    	gfArray.resize( sizeCounter );

    }




} /* namespace locust */


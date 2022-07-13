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
			fCavityFrequency( 1.0e9 ),
			fCavityQ( 1000 ),
			fThresholdFactor ( 0.25 ),
			fCavityDampingFactor( 0. ),
			fBFactor( 0. )
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

    bool DampedHarmonicOscillator::GenerateGreensFunction()
    {

        std::vector<std::pair<double,double>> tGFArray;

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
    	SetGFarray( tGFArray );


    	if ( tGFArray.size() < 1 ) return false;
    	else return true;

    }




} /* namespace locust */


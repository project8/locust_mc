/*
 * LMCGaussianNoiseGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGaussianNoiseGenerator.hh"

#include "logger.hh"
#include <random>


using std::string;

namespace locust
{
    LOGGER( lmclog, "GaussianNoiseGenerator" );

    MT_REGISTER_GENERATOR(GaussianNoiseGenerator, "gaussian-noise");

    GaussianNoiseGenerator::GaussianNoiseGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &GaussianNoiseGenerator::DoGenerateFreq ),
        fMean( 0. ),
        fSigma( 1. ),
        fRandomSeed(0),
        fNormDist( fMean, fSigma )
    {
        fRequiredSignalState = Signal::kFreq;
    }

    GaussianNoiseGenerator::~GaussianNoiseGenerator()
    {
    }

    bool GaussianNoiseGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "noise-floor" ) )
        {
            //tSigma = sqrt( aParam->get_value< double >( "noise-floor" ) * fAcquisitionRate * 1.e6);  // sampling rate
            fSigma = sqrt( aParam["noise-floor"]().as_double() );
            //printf("sigma is %g and fAcquisitionRate is %g\n", tSigma, fAcquisitionRate); getchar();
        }
        else
        {
            fSigma = aParam["sigma"]().as_double();
        }

  //      SetMeanAndSigma( aParam->get_value< double >( "mean", fMean ), tSigma );

        if (aParam.has( "random-seed") )
        {
            SetRandomSeed(  aParam.get_value< int >( "random-seed",fRandomSeed) );
        }


        if( aParam.has( "domain" ) )
        {
            string domain = aParam["domain"]().as_string();
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
                return false;
            }
        }

        return true;
    }

    void GaussianNoiseGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double GaussianNoiseGenerator::GetMean() const
    {
        return fMean;
    }

    void GaussianNoiseGenerator::SetMean( double aMean )
    {
        fMean = aMean;
        return;
    }

    double GaussianNoiseGenerator::GetSigma() const
    {
        return fSigma;
    }

    void GaussianNoiseGenerator::SetSigma( double aSigma )
    {
        fSigma = aSigma;
        return;
    }

    void GaussianNoiseGenerator::SetMeanAndSigma( double aMean, double aSigma, double aSampledSigma )
    {
        fNormDist = std::normal_distribution< double >( aMean, aSampledSigma );
        //fUniDist = std::uniform_real_distribution< double >(0.,360.);
        fMean = aMean;
        fSigma = aSigma;
        return;
    }
    int GaussianNoiseGenerator::GetRandomSeed() const
    {
        return fRandomSeed;
    }

    void GaussianNoiseGenerator::SetRandomSeed( int aRandomSeed )
    {
        fRandomSeed = aRandomSeed;
        return;
    }

    Signal::State GaussianNoiseGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void GaussianNoiseGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &GaussianNoiseGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &GaussianNoiseGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    bool GaussianNoiseGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool GaussianNoiseGenerator::DoGenerateTime( Signal* aSignal )
    {
        int random_seed_val;

        if ( fRandomSeed != 0 )
        {
            random_seed_val = fRandomSeed;
        }
        else
        {
            std::random_device rd;
            random_seed_val = rd();
        }
        std::default_random_engine generator(random_seed_val);

    	SetMeanAndSigma( fMean, fSigma, fSigma * sqrt(fAcquisitionRate * 1.e6) );

        double gain=1.;
        const unsigned nchannels = fNChannels;
        //double phi = 0.;  // voltage phase
        double mag_r = 0.;  // voltage mag
        double mag_i = 0.;

        for (int ch=0; ch<nchannels; ch++)
        {
            for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
            {
                //phi = fUniDist( fRNG );
                mag_r = fNormDist( generator ) * sqrt(0.5);
                mag_i = fNormDist( generator ) * sqrt(0.5);
                aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][0] += gain*sqrt(50.)* mag_r;
                aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][1] += gain*sqrt(50.)* mag_i;
                //	    printf("noise signal is %g\n", aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][1]); getchar();
            }
        }

        return true;
    }

    bool GaussianNoiseGenerator::DoGenerateFreq( Signal* aSignal )
    {
        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
        {
            aSignal->SignalFreq()[index][0] += fNormDist( fRNG );
            aSignal->SignalFreq()[index][1] += fNormDist( fRNG );
        }
        return true;
    }

} /* namespace locust */

/*
 * LMCGaussianNoiseGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGaussianNoiseGenerator.hh"

#include "logger.hh"

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
    fNormDist( fMean, fSigma )
    {
        fRequiredSignalState = Signal::kFreq;
    }

    GaussianNoiseGenerator::~GaussianNoiseGenerator()
    {
    }

    bool GaussianNoiseGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        double tSigma = fSigma;
        if( aParam->has( "noise-floor" ) )
        {
            //tSigma = sqrt( aParam->get_value< double >( "noise-floor" ) * aParam->get_value< double >( "acquisition-rate", 100. ) * 1.e6 );
            tSigma = sqrt( aParam->get_value< double >( "noise-floor" ) * 200.e6/2.);  // sampling rate is 200 MHz.  divide by 2 due to nbins->nbins/2 in katydid normalization.
        }
        else
        {
            tSigma = aParam->get_value( "sigma", fSigma );
        }

        SetMeanAndSigma( aParam->get_value< double >( "mean", fMean ), tSigma );

        if( aParam->has( "domain" ) )
        {
            string domain = aParam->get_value( "domain" );
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
        SetMeanAndSigma( aMean, fSigma );
        return;
    }

    double GaussianNoiseGenerator::GetSigma() const
    {
        return fSigma;
    }

    void GaussianNoiseGenerator::SetSigma( double aSigma )
    {
        SetMeanAndSigma( fMean, aSigma );
        return;
    }

    void GaussianNoiseGenerator::SetMeanAndSigma( double aMean, double aSigma )
    {
        fNormDist = std::normal_distribution< double >( aMean, aSigma );
        fMean = aMean;
        fSigma = aSigma;
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


    bool GaussianNoiseGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool GaussianNoiseGenerator::DoGenerateTime( Signal* aSignal ) const
    {
        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
            aSignal->SignalTime()[index] += fNormDist( fRNG );
        }

        return true;
    }

    bool GaussianNoiseGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
        {
            aSignal->SignalFreq()[index][0] += fNormDist( fRNG );
            aSignal->SignalFreq()[index][1] += fNormDist( fRNG );
        }
        return true;
    }

} /* namespace locust */

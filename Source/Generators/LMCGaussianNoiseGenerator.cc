/*
 * LMCGaussianNoiseGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGaussianNoiseGenerator.hh"

#include "Logger.hh"

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

    bool GaussianNoiseGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;

        double tSigma = fSigma;
        if( aParam->Has( "noise-floor" ) )
        {
            tSigma = sqrt( aParam->GetValue< double >( "noise-floor" ) * aParam->GetValue< double >( "acquisition-rate", 100. ) * 1.e6 );
        }
        else
        {
            tSigma = aParam->GetValue( "sigma", fSigma );
        }

        SetMeanAndSigma( aParam->GetValue< double >( "mean", fMean ), tSigma );

        if( aParam->Has( "domain" ) )
        {
            string domain = aParam->GetValue( "domain" );
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                DEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                ERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
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
            WARN( lmclog, "Unknown domain requested: " << aDomain );
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
            if (index<10) printf("signal %d is %g\n", index, aSignal->SignalTime()[index]);
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

/*
 * LMCGaussianNoiseGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGaussianNoiseGenerator.hh"

#include "../Core/LMCLogger.hh"

using std::string;

namespace locust
{
    LMCLOGGER( lmclog, "GaussianNoiseGenerator" );

    MT_REGISTER_GENERATOR(GaussianNoiseGenerator, "gaussian-noise");

    GaussianNoiseGenerator::GaussianNoiseGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &GaussianNoiseGenerator::DoGenerateFreq ),  // pls
            fMean( 0. ),
            fSigma( 1. ),
            fNormDist()
    {
        fRequiredSignalState = Signal::kFreq;  // pls
    }

    GaussianNoiseGenerator::~GaussianNoiseGenerator()
    {
    }

    bool GaussianNoiseGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;

        SetMean( aParam->GetValue< double >( "mean", fMean ) );
        SetSigma( aParam->GetValue< double >( "sigma", fSigma ) );

        if( aParam->Has( "domain" ) )
        {
            string domain = aParam->GetValue( "domain" );
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LMCDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LMCERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
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
            LMCWARN( lmclog, "Unknown domain requested: " << aDomain );
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
            aSignal->SignalTime( index ) += fNormDist( *fRNG, fMean, fSigma );
            if (index<10) printf("signal %d is %f\n", index, aSignal->SignalTime(index));
        }
        return true;
    }

    bool GaussianNoiseGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
        {
            aSignal->SignalFreq( index )[0] += fNormDist( *fRNG, fMean, fSigma );
            aSignal->SignalFreq( index )[1] += fNormDist( *fRNG, fMean, fSigma );
        }
        return true;
    }

} /* namespace locust */

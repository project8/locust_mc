/*
 * LMCGaussianNoiseGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGaussianNoiseGenerator.hh"

#include "LMCLogger.hh"

namespace locust
{
    LMCLOGGER( lmclog, "GaussianNoiseGenerator" );

    MT_REGISTER_GENERATOR(GaussianNoiseGenerator, "gaussian-noise");

    GaussianNoiseGenerator::GaussianNoiseGenerator( const std::string& aName ) :
            Generator( aName ),
            fMean( 0. ),
            fSigma( 1. ),
            fNormDist()
    {
        fRequiredSignalState = Signal::kTime;
    }

    GaussianNoiseGenerator::~GaussianNoiseGenerator()
    {
    }

    bool GaussianNoiseGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;

        SetMean( aParam->GetValue< double >( "mean", fMean ) );
        SetSigma( aParam->GetValue< double >( "sigma", fSigma ) );

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

    bool GaussianNoiseGenerator::DoGenerate( Signal* aSignal ) const
    {
        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
            aSignal->SignalTime( index ) += fNormDist( *fRNG, fMean, fSigma );
        }
        return true;
    }

} /* namespace locust */

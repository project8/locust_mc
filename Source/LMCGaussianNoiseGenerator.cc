/*
 * LMCGaussianNoiseGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGaussianNoiseGenerator.hh"

#include "LMCFactory.hh"

namespace locust
{
    static Registrar< Generator, GaussianNoiseGenerator > sGNGeneratorReg( "gaussian-noise" );

    GaussianNoiseGenerator::GaussianNoiseGenerator( const std::string& aName ) :
            Generator( aName ),
            fSigma( 1. ),
            fNormDist()
    {
        fRequiredSignalState = Signal::kTime;
    }

    GaussianNoiseGenerator::~GaussianNoiseGenerator()
    {
    }

    void GaussianNoiseGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return;

        SetSigma( aParam->get_value< double >( "sigma", fSigma ) );

        return;
    }

    void GaussianNoiseGenerator::Accept( GeneratorVisitor* aVisitor )
    {
        aVisitor->Visit( this );
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

    void GaussianNoiseGenerator::Generate( Signal* aSignal ) const
    {
        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
            aSignal->SignalTime( index ) += fNormDist( *fRNG, 0., fSigma );
        }
        return;
    }

} /* namespace locust */

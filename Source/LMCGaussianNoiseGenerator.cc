/*
 * LMCGaussianNoiseGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGaussianNoiseGenerator.hh"

namespace locust
{

    GaussianNoiseGenerator::GaussianNoiseGenerator() :
            Generator(),
            fSigma( 1. )
    {
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

        }
    }

} /* namespace locust */

/*
 * LMCGaussianNoiseGenerator.hh
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#ifndef LMCGAUSSIANNOISEGENERATOR_HH_
#define LMCGAUSSIANNOISEGENERATOR_HH_

#include "LMCGenerator.hh"
#include "LMCRunLengthCalculator.hh"

#include "RandomLib/NormalDistribution.hpp"

namespace locust
{

    /*!
     @class GaussianNoiseGenerator
     @author N. S. Oblath

     @brief Add Gaussian-distributed noise to the signal

     @details
     Can operate in time or frequency space

     Configuration name: "gaussian-noise"

     Available configuration options:
     - "noise-floor" -- Measured noise floor of system in W/Hz.  Overrides "sigma".
     - "mean": double -- Mean of the Gaussian noise (usually remains at 0)
     - "sigma": double -- Standard deviation of the Gaussian noise
     - "domain": string -- Determines whether the noise is generated in the time or frequency domain
                           Available options: "time" and "freq" [default]

    */
    class GaussianNoiseGenerator : public Generator
    {
        public:
            GaussianNoiseGenerator( const std::string& aName = "gaussian-noise" );
            virtual ~GaussianNoiseGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetNoiseFloor() const;
            void SetNoiseFloor( double aNoiseFloor );

            double GetMean() const;
            void SetMean( double aMean );

            double GetSigma() const;
            void SetSigma( double aSigma );

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );

        private:
            bool DoGenerate( Signal* aSignal ) const;

            bool DoGenerateTime( Signal* aSignal ) const;
            bool DoGenerateFreq( Signal* aSignal ) const;

            bool (GaussianNoiseGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;

            bool HasNoiseFloor;

            double fNoiseFloor;
            double fMean;
            double fSigma;
            RandomLib::NormalDistribution< double > fNormDist;
    };

} /* namespace locust */

#endif /* LMCGAUSSIANNOISEGENERATOR_HH_ */

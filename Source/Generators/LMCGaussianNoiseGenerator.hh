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
#include "LMCConst.hh"


#include <random>

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
     - "noise-floor" -- Measured noise floor of system in W/Hz.  Overrides "sigma"
     - "acquisition-rate" -- Digitization acquisition rate in MHz; Used in the calculation of sigma from the noise floor. Default is 100 MHz.
     - "mean": double -- Mean of the Gaussian noise (usually remains at 0)
     - "sigma": double -- Standard deviation of the Gaussian noise.  Ignored if "noise-floor" is present.
     - "domain": string -- Determines whether the noise is generated in the time or frequency domain
                           Available options: "time" and "freq" [default]

    */
    class GaussianNoiseGenerator : public Generator
    {
        public:
            GaussianNoiseGenerator( const std::string& aName = "gaussian-noise" );
            virtual ~GaussianNoiseGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetNoiseFloor() const;
            void SetNoiseFloor( double aNoiseFloor );

            double GetMean() const;
            void SetMean( double aMean );

            double GetSigma() const;
            void SetSigma( double aSigma );

            void SetMeanAndSigma( double aMean, double aSigma );

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );

        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (GaussianNoiseGenerator::*fDoGenerateFunc)( Signal* aSignal );

            double fMean;
            double fSigma;

            mutable std::normal_distribution< double > fNormDist;
            mutable std::uniform_real_distribution<double> fUniDist;


    };

} /* namespace locust */

#endif /* LMCGAUSSIANNOISEGENERATOR_HH_ */

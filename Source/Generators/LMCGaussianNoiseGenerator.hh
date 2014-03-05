/*
 * LMCGaussianNoiseGenerator.hh
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#ifndef LMCGAUSSIANNOISEGENERATOR_HH_
#define LMCGAUSSIANNOISEGENERATOR_HH_

#include "LMCGenerator.hh"

#include "RandomLib/NormalDistribution.hpp"

namespace locust
{

    /*!
     @class GaussianNoiseGenerator
     @author N. S. Oblath

     @brief Add Gaussian-distributed noise to the signal

     @details

     Configuration name: "gaussian-noise"

     Available configuration options:
     - "mean": double -- Mean of the Gaussian noise (usually remains at 0)
     - "sigma": double -- Standard deviation of the Gaussian noise

    */
    class GaussianNoiseGenerator : public Generator
    {
        public:
            GaussianNoiseGenerator( const std::string& aName = "gaussian-noise" );
            virtual ~GaussianNoiseGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetMean() const;
            void SetMean( double aMean );

            double GetSigma() const;
            void SetSigma( double aSigma );

        private:
            bool DoGenerate( Signal* aSignal ) const;

            double fMean;
            double fSigma;

            RandomLib::NormalDistribution< double > fNormDist;
    };

} /* namespace locust */

#endif /* LMCGAUSSIANNOISEGENERATOR_HH_ */

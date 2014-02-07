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

    class GaussianNoiseGenerator : public Generator
    {
        public:
            GaussianNoiseGenerator( const std::string& aName = "gaussian-noise" );
            virtual ~GaussianNoiseGenerator();

            void Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor );

            double GetSigma() const;
            void SetSigma( double aSigma );

        private:
            void Generate( Signal* aSignal ) const;

            double fSigma;

            RandomLib::NormalDistribution< double > fNormDist;
    };

} /* namespace locust */

#endif /* LMCGAUSSIANNOISEGENERATOR_HH_ */

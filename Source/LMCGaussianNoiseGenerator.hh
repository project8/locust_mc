/*
 * LMCGaussianNoiseGenerator.hh
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#ifndef LMCGAUSSIANNOISEGENERATOR_HH_
#define LMCGAUSSIANNOISEGENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    class GaussianNoiseGenerator : public Generator
    {
        public:
            GaussianNoiseGenerator();
            virtual ~GaussianNoiseGenerator();

            void Configure( const ParamNode* aNode );

            double GetSigma() const;
            void SetSigma( double aSigma );

        private:
            void Generate( Signal* aSignal ) const;

            double fSigma;
    };

} /* namespace locust */

#endif /* LMCGAUSSIANNOISEGENERATOR_HH_ */

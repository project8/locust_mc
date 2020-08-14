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

#ifdef ROOT_FOUND
    #include "LMCRunParameters.hh"
    #include "LMCRootTreeWriter.hh"
#endif

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
     - "noise-floor-psd": double -- Noise power in W/Hz.
     - "noise-temperature": double -- Noise temperature in K.
     - "domain": string -- Determines whether the noise is generated in the time or frequency domain
                           Available options: "time" and "freq" [default]
     - "write-root-tree": boolean -- Flag to control whether or not value of noise power is written
     	 	 			  to output Root file as a LMCRunParameter.
     - "random-seed": int -- Random seed used to generate random noise.  If this is omitted then the
     	 	 	 	 	  noise spectrum is reproducible.
     - "root-filename": string -- Name of output Root file.  This can have the same name as other
          	 	 	 	  generators' output Root files, in which case all of the Root objects will
          	 	 	 	  be written to the same output file.

    */
    class GaussianNoiseGenerator : public Generator
    {
        public:
            GaussianNoiseGenerator( const std::string& aName = "gaussian-noise" );
            virtual ~GaussianNoiseGenerator();

            bool Configure( const scarab::param_node& aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetNoiseFloor() const;
            void SetNoiseFloor( double aNoiseFloor );

            double GetMean() const;
            void SetMean( double aMean );

            double GetSigma() const;
            void SetSigma( double aSigma );

            void SetMeanAndSigma( double aMean, double aSigma, double aSampledSigma);

            int GetRandomSeed() const;
            void SetRandomSeed(  int aRandomSeed );


            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );

        private:
            bool WriteRootTree();
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (GaussianNoiseGenerator::*fDoGenerateFunc)( Signal* aSignal );

            double fMean;
            double fSigma;
            int fRandomSeed;
            bool fWriteRootTree;
            std::string fRootFilename;

            mutable std::normal_distribution< double > fNormDist;
    };

} /* namespace locust */

#endif /* LMCGAUSSIANNOISEGENERATOR_HH_ */

/*
 * LMCHighPassFilterFFTGenerator.hh
 *
 *  Created on: 14 Feb 2019
 *      Author: plslocum
 */

#ifndef LMCHIGHPASSFILTERFFTGENERATOR_HH_
#define LMCHIGHPASSFILTERFFTGENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class HighPassFilterFFTGenerator
     @author P. L. Slocum

     @brief Apply high pass filter to signal in frequency domain

     @details
      This is a filter in frequency space.  It assumes that the signal has been
      sampled at fAcquisitionRate*aSignal->DecimationFactor(), and needs to be
      followed by the decimation generator.
      Available configuration options:
      - "threshold": double -- threshold below which signals are suppressed (Hz).

     Configuration name: "hpf-fft"

    */

class HighPassFilterFFTGenerator : public Generator
    {
        public:
            HighPassFilterFFTGenerator( const std::string& aName = "hpf-fft" );
            virtual ~HighPassFilterFFTGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;
            double GetThreshold() const;
            void SetThreshold( double aThreshold );

        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (HighPassFilterFFTGenerator::*fDoGenerateFunc)( Signal* aSignal );
            double fThreshold; // (Hz)

    };

} /* namespace locust */

#endif /* LMCLOWPASSFILTERFFTGENERATOR_HH_ */

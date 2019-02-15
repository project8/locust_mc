/*
 * LMCLowPassFilterFFTGenerator.hh
 *
 *  Created on: 29 January 2015
 *      Author: plslocum after nsoblath
 */

#ifndef LMCLOWPASSFILTERFFTGENERATOR_HH_
#define LMCLOWPASSFILTERFFTGENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class LowPassFilterFFTGenerator
     @author P. L. Slocum after N. S. Oblath

     @brief Apply low pass filter to fast-sampled signal using an fft.

     @details
     This filter needs to be followed by decimation.  It assumes that the signal has been sampled faster (X aSignal->DecimationFactor() ) than the desired sampling frequency so that high frequency spurs are measured correctly.

      Available configuration options:
      - "threshold": double -- ratio of cutoff frequency to Nyquist frequency


     Configuration name: "lpf-fft"

    */

class LowPassFilterFFTGenerator : public Generator
    {
        public:
            LowPassFilterFFTGenerator( const std::string& aName = "lpf-fft" );
            virtual ~LowPassFilterFFTGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;
            double GetThreshold() const;
            void SetThreshold( double aThreshold );



        private:
            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );
            bool DoGenerateFreq( Signal* aSignal );

            bool (LowPassFilterFFTGenerator::*fDoGenerateFunc)( Signal* aSignal );

            double fThreshold;

    };

} /* namespace locust */

#endif /* LMCLOWPASSFILTERFFTGENERATOR_HH_ */

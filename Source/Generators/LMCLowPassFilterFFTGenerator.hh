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
     Can operate in frequency or time space, but only the frequency version will have any effect on the signal.

     Configuration name: "lpf-fft"

    */

class LowPassFilterFFTGenerator : public Generator
    {
        public:
            LowPassFilterFFTGenerator( const std::string& aName = "lpf-fft" );
            virtual ~LowPassFilterFFTGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetReceiverGain() const;
            void SetReceiverGain( double aReceiverGain );


        private:
            bool DoGenerate( Signal* aSignal ) const;

            bool DoGenerateTime( Signal* aSignal ) const;
            bool DoGenerateFreq( Signal* aSignal ) const;

            bool (LowPassFilterFFTGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;

            double fReceiverGain;

    };

} /* namespace locust */

#endif /* LMCLOWPASSFILTERFFTGENERATOR_HH_ */

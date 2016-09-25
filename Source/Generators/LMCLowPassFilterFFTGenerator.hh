/*
 * LMCLowPassFilterGenerator.hh
 *
 *  Created on: 29 January 2015
 *      Author: plslocum after nsoblath
 */

#ifndef LMCLOWPASSFILTERGENERATOR_HH_
#define LMCLOWPASSFILTERGENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class LowPassFilterGenerator
     @author P. L. Slocum after N. S. Oblath

     @brief Apply low pass filter to signal using an fft.

     @details
     Can operate in frequency or time space, but only the frequency version will have any effect on the signal.

     Configuration name: "lpf-fft"

    */

class LowPassFilterGenerator : public Generator
    {
        public:
            LowPassFilterGenerator( const std::string& aName = "lpf-fft" );
            virtual ~LowPassFilterGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetReceiverGain() const;
            void SetReceiverGain( double aReceiverGain );


        private:
            bool DoGenerate( Signal* aSignal ) const;

            bool DoGenerateTime( Signal* aSignal ) const;
            bool DoGenerateFreq( Signal* aSignal ) const;

            bool (LowPassFilterGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;

            double fReceiverGain;

    };

} /* namespace locust */

#endif /* LMCLOWPASSFILTERGENERATOR_HH_ */

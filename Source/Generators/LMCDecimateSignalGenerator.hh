/*
 * LMCDecimateSignalGenerator.hh
 *
 *  Created on: 21 September 2016
 *      Author: plslocum after nsoblath
 */

#ifndef LMCDECIMATESIGNALGENERATOR_HH_
#define LMCDECIMATESIGNALGENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class DecimateSignalGenerator
     @author P. L. Slocum after N. S. Oblath

     @brief Decimate signal to lower the effective sampling frequency.

     @details
     Operates in the time domain.

     Configuration name: "decimate"

    */

class DecimateSignalGenerator : public Generator
    {
        public:
            DecimateSignalGenerator( const std::string& aName = "decimate-signal" );
            virtual ~DecimateSignalGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetReceiverGain() const;
            void SetReceiverGain( double aReceiverGain );


        private:
            bool DoGenerate( Signal* aSignal ) const;

            bool DoGenerateTime( Signal* aSignal ) const;
            bool DoGenerateFreq( Signal* aSignal ) const;

            bool (DecimateSignalGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;

            double fReceiverGain;

    };

} /* namespace locust */

#endif /* LMCDECIMATESIGNALGENERATOR_HH_ */

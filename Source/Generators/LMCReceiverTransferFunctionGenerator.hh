/*
 * LMCReceiverTransferFunctionGenerator.hh
 *
 *  Created on: 29 January 2015
 *      Author: plslocum after nsoblath
 */

#ifndef LMCRECEIVERTRANSFERFUNCTIONGENERATOR_HH_
#define LMCRECEIVERTRANSFERFUNCTIONGENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class ReceiverTransferFunctionGenerator
     @author P. L. Slocum after N. S. Oblath

     @brief Apply receiver transfer function to signal.

     @details
     Can operate in frequency or time space, but only the frequency version will have any effect on the signal.

     Configuration name: "receiver-transfer-function"

     Available configuration options:
     - "receiver-gain" double -- total gain of receiver at one frequency after downmixing (dB).

    */
    class ReceiverTransferFunctionGenerator : public Generator
    {
        public:
            ReceiverTransferFunctionGenerator( const std::string& aName = "receiver-transfer-function" );
            virtual ~ReceiverTransferFunctionGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetReceiverGain() const;
            void SetReceiverGain( double aReceiverGain );


        private:
            bool DoGenerate( Signal* aSignal ) const;

            bool DoGenerateTime( Signal* aSignal ) const;
            bool DoGenerateFreq( Signal* aSignal ) const;

            bool (ReceiverTransferFunctionGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;

            double fReceiverGain;

    };

} /* namespace locust */

#endif /* LMCRECEIVERTRANSFERFUNCTIONGENERATOR_HH_ */

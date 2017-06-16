/*
 * LMCButterworthLPFGenerator.hh
 *
 *  Created on: 29 January 2015
 *      Author: plslocum after nsoblath
 */

#ifndef LMCBUTTERWORTHLPFGENERATOR_HH_
#define LMCBUTTERWORTHLPFGENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class ButterworthLPFGenerator
     @author P. L. Slocum after N. S. Oblath

     @brief Apply digital low pass filter to signal.

     @details
     Can operate in frequency or time space, but only the frequency version will have any effect on the signal.

     Configuration name: "butterworth-lpf"

    */

class ButterworthLPFGenerator : public Generator
    {
        public:
            ButterworthLPFGenerator( const std::string& aName = "butterworth-lpf" );
            virtual ~ButterworthLPFGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetReceiverGain() const;
            void SetReceiverGain( double aReceiverGain );


        private:
            bool DoGenerate( Signal* aSignal ) const;

            bool DoGenerateTime( Signal* aSignal ) const;
            bool DoGenerateFreq( Signal* aSignal ) const;

            bool (ButterworthLPFGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;

            double fReceiverGain;

    };

} /* namespace locust */

#endif /* LMCBUTTERWORTHLPFGENERATOR_HH_ */

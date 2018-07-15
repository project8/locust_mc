/*
 * LMCLocalOscillatorGenerator.hh
 *
 *  Created on: Apr 20 2018
 *      Author: buzinsky
 */

#ifndef LMCLOCALOSCILLATORGENERATOR_HH_
#define LMCLOCALOSCILLATORGENERATOR_HH_

#include "LMCGenerator.hh"

namespace locust
{

    /*!
     @class LocalOscillatorGenerator
     @author N. Buzinsky

     @brief Downmix the signal (multiplicative) by f_LO

     @details
     Time domain only

     Configuration name: "local-oscillator"

    */

class LocalOscillatorGenerator : public Generator
    {
        public:
            LocalOscillatorGenerator( const std::string& aName = "local-oscillator" );
            virtual ~LocalOscillatorGenerator();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

        private:
            double fLOFrequency;
            double* ComplexMultiplication(fftw_complex a, fftw_complex b);

            bool DoGenerate( Signal* aSignal );

            bool DoGenerateTime( Signal* aSignal );

            bool (LocalOscillatorGenerator::*fDoGenerateFunc)( Signal* aSignal );
    };

} /* namespace locust */

#endif /* LMCLOCALOSCILLATORGENERATOR_HH_ */

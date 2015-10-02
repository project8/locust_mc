/*
 * LMCKassiopeiaGenerator.hh
 *
 *  Created on: Feb 4, 2014
 *      Author: plslocum after nsoblath
 */

#ifndef LMCKASSIOPEIAGENERATOR_HH_
#define LMCKASSIOPEIAGENERATOR_HH_

#include "LMCGenerator.hh"
#include "LMCRunLengthCalculator.hh"
#include <pthread.h>


namespace locust
{

    /*!
     @class KassiopeiaGenerator
     @author P. L. Slocum after N. S. Oblath

     @brief Interface with Kassiopeia

    */
    class KassiopeiaGenerator : public Generator
    {
        public:
            KassiopeiaGenerator( const std::string& aName = "kassiopeia" );
            virtual ~KassiopeiaGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );

        private:
            bool DoGenerate( Signal* aSignal ) const;

//            void * KassiopeiaInit(void *);

            bool DoGenerateTime( Signal* aSignal ) const;
            bool DoGenerateFreq( Signal* aSignal ) const;

            bool (KassiopeiaGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;
    };

} /* namespace locust */

#endif /* LMCKASSIOPEIAGENERATOR_HH_ */

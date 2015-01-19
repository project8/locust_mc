/*
 * LMCTestSignalGenerator.hh
 *
 *  Created on: Dec 15, 2014
 *      Author: pslocum
 */

#ifndef LMCTESTSIGNALGENERATOR_HH_
#define LMCTESTSIGNALGENERATOR_HH_

#include "../Core/LMCGenerator.hh"

namespace locust
{

    /*!
     @class [name]Generator
     @author N. S. Oblath

     @brief

     @details
     Operates in [domain] space

     Configuration name: "config-name"

     Available configuration options:
     - "param-name": type -- Description

    */
    class TestSignalGenerator : public Generator
    {
        public:
            TestSignalGenerator( const std::string& aName = "test-signal" );
            virtual ~TestSignalGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

        private:
            bool DoGenerate( Signal* aSignal ) const;

            bool DoGenerateTime( Signal* aSignal ) const;
//            bool DoGenerateFreq( Signal* aSignal ) const;

            bool (TestSignalGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;


    };

} /* namespace locust */

#endif /* LMCTESTSIGNALGENERATOR_HH_ */

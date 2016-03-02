/*
 * LMCGenerator.hh
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#ifndef LMCGENERATOR_HH_
#define LMCGENERATOR_HH_

#include "LMCFactory.hh"
#include "LMCVisitor.hh"
#include "LMCParam.hh"
#include "LMCSignal.hh"

#include <random>

namespace locust
{

    class Generator
    {
        public:
            Generator( const std::string& aName = "generic-generator" );
            virtual ~Generator();

            virtual bool Configure( const ParamNode* aNode ) = 0;

            virtual void Accept( GeneratorVisitor* aVisitor ) const = 0;

            Signal* Run( unsigned aTimeSize ) const;
            bool Run( Signal* aSignal ) const;

            const std::string& GetName() const;
            void SetName( const std::string& aName );

            Signal::State GetRequiredSignalState() const;
            void SetRequiredSignalState( Signal::State state );

            Generator* GetNextGenerator() const;
            void SetNextGenerator( Generator* aGenerator );

            static std::mt19937_64& RNG();

        protected:
            bool Generate( Signal* aSignal ) const;
            virtual bool DoGenerate( Signal* aSignal ) const = 0;

            std::string fName;

            Signal::State fRequiredSignalState;

            Generator* fNext;

            static std::mt19937_64 fRNG;
    };


#define MT_REGISTER_GENERATOR(gen_class, gen_name) \
        static Registrar< Generator, gen_class > s_##gen_class##_writer_registrar( gen_name );

} /* namespace locust */

#endif /* LMCGENERATOR_HH_ */

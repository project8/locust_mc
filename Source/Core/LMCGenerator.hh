/*
 * LMCGenerator.hh
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#ifndef LMCGENERATOR_HH_
#define LMCGENERATOR_HH_

#include "LMCVisitor.hh"
#include "param.hh"
#include "LMCSignal.hh"

#include "factory.hh"

#include <random>

namespace locust
{

    class Generator
    {
        public:
            Generator( const std::string& aName = "generic-generator" );
            virtual ~Generator();

            virtual bool Configure( const scarab::param_node& aNode ) = 0;
            void ConfigureNChannels( unsigned nchannels ) const;
            void ConfigureAcquisitionRate( double ar ) const;

            virtual void Accept( GeneratorVisitor* aVisitor ) const = 0;

            Signal* Run( unsigned aTimeSize );
            bool Run( Signal* aSignal );

            const std::string& GetName() const;
            void SetName( const std::string& aName );

            Signal::State GetRequiredSignalState() const;
            void SetRequiredSignalState( Signal::State state );

            Generator* GetNextGenerator() const;
            void SetNextGenerator( Generator* aGenerator );

      //            bool PrepareNChannels( const Digitizer* aDig ) const;

            static std::mt19937_64& RNG();

        protected:
            bool Generate( Signal* aSignal );
            virtual bool DoGenerate( Signal* aSignal ) = 0;
            mutable unsigned fNChannels;
            mutable double fAcquisitionRate;


            std::string fName;

            Signal::State fRequiredSignalState;

            Generator* fNext;

            static std::mt19937_64 fRNG;
    };


#define MT_REGISTER_GENERATOR(gen_class, gen_name) \
        static scarab::registrar< Generator, gen_class > s_##gen_class##_generator_registrar( gen_name );

} /* namespace locust */

#endif /* LMCGENERATOR_HH_ */

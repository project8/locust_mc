/*
 * LMCGenerator.hh
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#ifndef LMCGENERATOR_HH_
#define LMCGENERATOR_HH_

#include "LMCVisitor.hh"
#include "LMCParam.hh"
#include "LMCSignal.hh"

#include "RandomLib/Random.hpp"

namespace locust
{

    class Generator
    {
        public:
            Generator( const std::string& aName = "generic-generator" );
            virtual ~Generator();

            virtual void Configure( const ParamNode* aNode ) = 0;

            virtual void Accept( GeneratorVisitor* aVisitor ) const = 0;

            Signal* Run( unsigned aTimeSize ) const;
            void Run( Signal* aSignal ) const;

            const std::string& GetName() const;
            void SetName( const std::string& aName );

            Signal::State GetRequiredSignalState() const;
            void SetRequiredSignalState( Signal::State state );

            void SetRNG( RandomLib::Random* aRNG );
            RandomLib::Random* GetRNG() const;

            Generator* GetNextGenerator() const;
            void SetNextGenerator( Generator* aGenerator );

        protected:
            virtual void Generate( Signal* aSignal ) const = 0;

            std::string fName;

            Signal::State fRequiredSignalState;

            RandomLib::Random* fRNG;

            Generator* fNext;
    };

} /* namespace locust */

#endif /* LMCGENERATOR_HH_ */

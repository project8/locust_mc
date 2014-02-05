/*
 * LMCGenerator.hh
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#ifndef LMCGENERATOR_HH_
#define LMCGENERATOR_HH_

#include "LMCParam.hh"
#include "LMCSignal.hh"

#include "RandomLib/Random.hpp"

namespace locust
{

    class Generator
    {
        public:
            Generator();
            virtual ~Generator();

            virtual void Configure( const ParamNode* aNode ) = 0;

            void Run( unsigned aTimeSize ) const;
            void Run( Signal* aSignal ) const;

            Signal::State GetRequiredSignalState() const;
            void SetRequiredSignalState( Signal::State state );

            void SetRNG( RandomLib::Random* aRNG );
            RandomLib::Random* GetRNG() const;

            const Generator* GetNextGenerator() const;
            void SetNextGenerator( const Generator* aGenerator );

        protected:
            virtual void Generate( Signal* aSignal ) const = 0;

            Signal::State fRequiredSignalState;

            RandomLib::Random* fRNG;

            const Generator* fNext;
    };

} /* namespace locust */

#endif /* LMCGENERATOR_HH_ */

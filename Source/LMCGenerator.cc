/*
 * LMCGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGenerator.hh"

#include "LMCSignal.hh"

namespace locust
{

    Generator::Generator() :
            fRequiredSignalState( Signal::kTime ),
            fRNG( NULL ),
            fNext( NULL )
    {
    }

    Generator::~Generator()
    {
    }

    void Generator::Run( unsigned aTimeSize ) const
    {
        Signal* newSignal = new Signal();
        newSignal->Initialize( aTimeSize );
        Run( newSignal );
        return;
    }

    void Generator::Run( Signal* aSignal ) const
    {
        aSignal->ToState( fRequiredSignalState );
        Generate( aSignal );
        if( fNext != NULL ) fNext->Run( aSignal );
        return;
    }

    Signal::State Generator::GetRequiredSignalState() const
    {
        return fRequiredSignalState;
    }

    void Generator::SetRequiredSignalState( Signal::State state )
    {
        fRequiredSignalState = state;
        return;
    }

    void Generator::SetRNG( RandomLib::Random* aRNG )
    {
        fRNG = aRNG;
        return;
    }

    RandomLib::Random* Generator::GetRNG() const
    {
        return fRNG;
    }

    const Generator* Generator::GetNextGenerator() const
    {
        return fNext;
    }

    void Generator::SetNextGenerator( const Generator* aGenerator )
    {
        fNext = aGenerator;
        return;
    }


} /* namespace locust */

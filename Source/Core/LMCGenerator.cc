/*
 * LMCGenerator.cc
 *
 *  Created on: Feb 4, 2014
 *      Author: nsoblath
 */

#include "LMCGenerator.hh"

#include "logger.hh"
#include "LMCSignal.hh"

namespace locust
{
    LOGGER( lmclog, "Generator" );

    std::mt19937_64 Generator::fRNG;

    std::mt19937_64& Generator::RNG()
    {
        return Generator::fRNG;
    }

    Generator::Generator( const std::string& aName ) :
            fName( aName ),
            fRequiredSignalState( Signal::kFreq ),
            fNChannels( 1 ),
            fNext( NULL )
    {
    }

    Generator::~Generator()
    {
    }

    Signal* Generator::Run( unsigned aTimeSize ) 
    {
        Signal* newSignal = new Signal();
        newSignal->Initialize( aTimeSize , fNChannels );
        if( ! Run( newSignal ) )
        {
            delete newSignal;
            return NULL;
        }
        return newSignal;
    }

    bool Generator::Run( Signal* aSignal )
    {
        if(! Generate( aSignal ) )
        {
            LERROR( lmclog, "Signal generation failed" );
            return false;
        }
        if( fNext != NULL ) fNext->Run( aSignal );
        return true;
    }

    bool Generator::Generate( Signal* aSignal )
    {
        if( ! aSignal->ToState( fRequiredSignalState ) )
        {
            LERROR( lmclog, "Unable to convert signal to state <" << fRequiredSignalState << ">" );
            return false;
        }
        return DoGenerate( aSignal );
    }

    const std::string& Generator::GetName() const
    {
        return fName;
    }

    void Generator::SetName( const std::string& aName )
    {
        fName = aName;
        return;
    }

  void Generator::ConfigureNChannels( unsigned nchannels ) const
  {
    fNChannels = nchannels;
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

    Generator* Generator::GetNextGenerator() const
    {
        return fNext;
    }

    void Generator::SetNextGenerator( Generator* aGenerator )
    {
        fNext = aGenerator;
        return;
    }


} /* namespace locust */

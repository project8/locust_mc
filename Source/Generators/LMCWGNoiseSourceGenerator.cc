/*
 * LMCWGNoiseSourceGenerator.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#include "LMCWGNoiseSourceGenerator.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "WGNoiseSourceGenerator" );

    MT_REGISTER_GENERATOR(WGNoiseSourceGenerator, "noise-source");

    WGNoiseSourceGenerator::WGNoiseSourceGenerator( const std::string& aName ) :
            Generator( aName )
    {
        fRequiredSignalState = Signal::kFreq;
    }

    WGNoiseSourceGenerator::~WGNoiseSourceGenerator()
    {
    }

    bool WGNoiseSourceGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        return true;
    }

    void WGNoiseSourceGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    bool WGNoiseSourceGenerator::DoGenerate( Signal* aSignal ) const
    {
        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
        {
            aSignal->SignalFreq()[index] += ???;
        }
        return true;
    }

    void WGNoiseSourceGenerator::GenerateWGe( unsigned aFreqSize, double aFreqBW )
    {
        fWGe.clear();
        fWGe.resize( aFreqSize );

    }

} /* namespace locust */

/*
 * LMCDecimateSignalGenerator.cc
 *
 *  Created on: Sept 9, 2016
 *      Author: plslocum after nsoblath
 */

#include "LMCDecimateSignalGenerator.hh"

#include "logger.hh"
#include "LMCGlobalsDeclaration.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "DecimateSignalGenerator" );

    MT_REGISTER_GENERATOR(DecimateSignalGenerator, "decimate-signal");

    DecimateSignalGenerator::DecimateSignalGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &DecimateSignalGenerator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::kTime;
    }

    DecimateSignalGenerator::~DecimateSignalGenerator()
    {
    }

    bool DecimateSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;
        return true;
    }

    void DecimateSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    bool DecimateSignalGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool DecimateSignalGenerator::DoGenerateTime( Signal* aSignal ) const
    {
        // Decimate Fs -> Fs/10
        const int fDecimationFactor = 10;
        for (int ch=0; ch<NCHANNELS; ch++)
        {
            for( unsigned index = 0; index < aSignal->TimeSize()*fDecimationFactor; ++index )
            {
                if (index % fDecimationFactor == 0)
                {
                    aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index/fDecimationFactor][0] =
                    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*fDecimationFactor + index][0];
                    aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index/fDecimationFactor][1] =
                    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*fDecimationFactor + index][1];
                }
            }
        }

         return true;
    }

    bool DecimateSignalGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
        return true;
    }

} /* namespace locust */

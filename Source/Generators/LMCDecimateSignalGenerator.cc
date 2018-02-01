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


    bool DecimateSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool DecimateSignalGenerator::DoGenerateTime( Signal* aSignal ) 
    {
        // Decimate Fs -> Fs/DecimationFactor
        for (int ch=0; ch<NCHANNELS; ++ch)
        {
            for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
            {
                if (index % aSignal->DecimationFactor() == 0)
                {
                    aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index/aSignal->DecimationFactor()][0] =
                    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0];
                    aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index/aSignal->DecimationFactor()][1] =
                    aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1];
                }
            }
        }

         return true;
    }

} /* namespace locust */

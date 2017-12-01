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
/*


    double* aTemporarySignal = new double[aSignal->TimeSize()];

    // first copy then zero the Signal.
    for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
      {
      aTemporarySignal[index] = aSignal->SignalTime()[index];
      aSignal->SignalTime()[index] = 0.;
      }

      */

    // Decimate Fs -> Fs/10
    for( unsigned index = 0; index < 41943040; ++index )
      {
      if (index%fDecimationFactor == 0)
        {
        aSignal->SignalTime()[index/fDecimationFactor] = aLongSignal[index];

        aSignal->SignalTimeComplex()[index/fDecimationFactor][0] =
        		aSignal->LongSignalTimeComplex()[index][0];
        aSignal->SignalTimeComplex()[index/fDecimationFactor][1] =
        		aSignal->LongSignalTimeComplex()[index][1];
        }
      }

//    delete aTemporarySignal;

         return true;
    }

    bool DecimateSignalGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
        return true;
    }

} /* namespace locust */

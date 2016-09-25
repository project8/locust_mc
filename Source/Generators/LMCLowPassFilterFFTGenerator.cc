/*
 * LMCLowPassFilterGenerator.cc
 *
 *  Created on: Sept 9, 2016
 *      Author: plslocum after nsoblath
 */

#include "LMCLowPassFilterFFTGenerator.hh"

#include "logger.hh"

using std::string;

namespace locust
{
    LOGGER( lmclog, "LowPassFilterGenerator" );

    MT_REGISTER_GENERATOR(LowPassFilterGenerator, "lpf-fft");

    LowPassFilterGenerator::LowPassFilterGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &LowPassFilterGenerator::DoGenerateFreq )
    {
        fRequiredSignalState = Signal::kFreq;
    }

    LowPassFilterGenerator::~LowPassFilterGenerator()
    {
    }

    bool LowPassFilterGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;
        return true;
    }

    void LowPassFilterGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    bool LowPassFilterGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool LowPassFilterGenerator::DoGenerateTime( Signal* aSignal ) const
    {
    return true;
    }

    bool LowPassFilterGenerator::DoGenerateFreq( Signal* aSignal ) const
    {


      double CutoffFreq = 85.e6;

        for( unsigned index = 0; index < aSignal->FreqSize(); ++index )
          {
	    //	    printf("index is %d and cutofffreq is %g\n", index, aSignal->FreqSize()*CutoffFreq/1.e9);
        	// LPF
        	if (index > aSignal->FreqSize()*CutoffFreq/1.e9)
        	{
          	  aSignal->SignalFreq()[index][0] = 0.;
        	  aSignal->SignalFreq()[index][1] = 0.;
        	}
        	// normalize
        	aSignal->SignalFreq()[index][0] = aSignal->SignalFreq()[index][0]/(double)aSignal->TimeSize();
        	aSignal->SignalFreq()[index][1] = aSignal->SignalFreq()[index][1]/(double)aSignal->TimeSize();

          }



        return true;
    }

} /* namespace locust */

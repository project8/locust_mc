/*
 * LMCLowPassFilterGenerator.cc
 *
 *  Created on: Sept 9, 2016
 *      Author: plslocum after nsoblath
 */

#include "LMCLowPassFilterFFTGenerator.hh"

#include "logger.hh"
#include "LMCGlobalsDeclaration.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "LowPassFilterGenerator" );

    MT_REGISTER_GENERATOR(LowPassFilterGenerator, "lpf-fft");

    LowPassFilterGenerator::LowPassFilterGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &LowPassFilterGenerator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::kTime;
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


        double CutoffFreq = 85.e6;
        int nwindows = 80;
        int windowsize = 10*aSignal->TimeSize()/nwindows;


      	fftw_complex *SignalComplex;
          SignalComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );
      	fftw_complex *FFTComplex;
          FFTComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );

          fftw_plan ForwardPlan;
          ForwardPlan = fftw_plan_dft_1d(windowsize, SignalComplex, FFTComplex, FFTW_FORWARD, FFTW_ESTIMATE);
          fftw_plan ReversePlan;
          ReversePlan = fftw_plan_dft_1d(windowsize, FFTComplex, SignalComplex, FFTW_BACKWARD, FFTW_ESTIMATE);


          for (int nwin = 0; nwin < nwindows; nwin++)
          {


      // Construct complex voltage.
          for( unsigned index = 0; index < windowsize; ++index )
          {
              SignalComplex[index][0] = aLongSignal[ nwin*windowsize + index ];
              SignalComplex[index][1] = 0.;
      //        if (index==20000) {printf("signal 20000 is %g\n", aSignal->SignalTime()[index]); getchar();}

          }

          fftw_execute(ForwardPlan);


// Low Pass Filter
        for( unsigned index = 0; index < windowsize; ++index )
          {
          if
            (
          	(index > windowsize/2.*CutoffFreq/1.e9)
            )
        	{
//        	  printf("index is %d\n", index);
              FFTComplex[index][0] = 0.;
              FFTComplex[index][1] = 0.;
            }
          }


          fftw_execute(ReversePlan);

          double norm = (double)(windowsize);

          for( unsigned index = 0; index < windowsize; ++index )
          {
          	// normalize and take the real part of the reverse transform, for digitization.
      //      aSignal->SignalTime()[ nwin*windowsize + index ] = SignalComplex[index][0]/norm;
            aLongSignal[ nwin*windowsize + index ] = SignalComplex[index][0]/norm;
            //    if (index>=20000) {printf("filtered signal is %g\n", aSignal->SignalTime()[index]); getchar();}
          }


          }  // nwin

      delete SignalComplex;
      delete FFTComplex;



    	return true;
    }

    bool LowPassFilterGenerator::DoGenerateFreq( Signal* aSignal ) const
    {


        return true;
    }

} /* namespace locust */

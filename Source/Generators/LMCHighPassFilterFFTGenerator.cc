/*
 * LMCHighPassFilterFFTGenerator.cc
 *
 *  Created on: Feb. 14 2019
 *      Author: plslocum
 */

#include "LMCHighPassFilterFFTGenerator.hh"

#include "logger.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "HighPassFilterFFTGenerator" );

    MT_REGISTER_GENERATOR(HighPassFilterFFTGenerator, "hpf-fft");

    HighPassFilterFFTGenerator::HighPassFilterFFTGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &HighPassFilterFFTGenerator::DoGenerateTime ),
        fThreshold ( 2.0e6 ) // MHz
    {
        fRequiredSignalState = Signal::kTime;
    }

    HighPassFilterFFTGenerator::~HighPassFilterFFTGenerator()
    {
    }

    bool HighPassFilterFFTGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "threshold" ) )
        {
            SetThreshold( aParam.get_value< double >( "threshold", fThreshold ) ); 
        }
        return true;
    }

    double HighPassFilterFFTGenerator::GetThreshold() const
    {
        return fThreshold;
    }

    void HighPassFilterFFTGenerator::SetThreshold( double aThreshold )
    {
        fThreshold = aThreshold;
        return;
    }


    void HighPassFilterFFTGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    bool HighPassFilterFFTGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool HighPassFilterFFTGenerator::DoGenerateTime( Signal* aSignal )
    {

        const unsigned nchannels = fNChannels;


        double CutoffFreq = fThreshold; // Hz
        double FastNyquist = fAcquisitionRate/2. * 1.e6 * aSignal->DecimationFactor();
        int nwindows = 80;
        int windowsize = aSignal->DecimationFactor()*aSignal->TimeSize()/nwindows;


        fftw_complex *SignalComplex;
        SignalComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );
        fftw_complex *FFTComplex;
        FFTComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );

        fftw_plan ForwardPlan;
        ForwardPlan = fftw_plan_dft_1d(windowsize, SignalComplex, FFTComplex, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan ReversePlan;
        ReversePlan = fftw_plan_dft_1d(windowsize, FFTComplex, SignalComplex, FFTW_BACKWARD, FFTW_ESTIMATE);

        for (int ch=0; ch<nchannels; ch++)
        {
            for (int nwin = 0; nwin < nwindows; nwin++)
            {
                // Construct complex voltage.
                for( unsigned index = 0; index < windowsize; ++index )
                {
                    SignalComplex[index][0] = aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + nwin*windowsize + index ][0];
                    SignalComplex[index][1] = aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + nwin*windowsize + index ][1];
                    //if (index==20000) {printf("signal 20000 is %g\n", aSignal->SignalTime()[index]); getchar();}
                }

                fftw_execute(ForwardPlan);

//                printf("windowsize/2 is %g and CutoffFreq is %g and FastNyquist is %g\n", windowsize/2, CutoffFreq, FastNyquist); getchar();

                // High Pass FilterFFT
                for( unsigned index = 0; index < windowsize; ++index )
                {
                if ((index < windowsize/2. * CutoffFreq/FastNyquist)||(index > windowsize/2.))
                  {
                  FFTComplex[index][0] = 0.;
                  FFTComplex[index][1] = 0.;
                  }
                }

                fftw_execute(ReversePlan);

                double norm = (double)(windowsize);

                for( unsigned index = 0; index < windowsize; ++index )
                {
                    // normalize
                    aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + nwin*windowsize + index ][0] = SignalComplex[index][0]/norm;
                    aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + nwin*windowsize + index ][1] = SignalComplex[index][1]/norm;
                    //if (index>=20000) {printf("filtered signal is %g\n", aSignal->SignalTime()[index]); getchar();}
                }
            }  // nwin
        }  // NCHANNELS

        delete [] SignalComplex;
        delete [] FFTComplex;

        return true;
    }

    bool HighPassFilterFFTGenerator::DoGenerateFreq( Signal* aSignal )
    {
       return true;
    }

} /* namespace locust */

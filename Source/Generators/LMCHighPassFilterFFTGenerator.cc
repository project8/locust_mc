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
        fThreshold ( 1.0e6 ) // 1 MHz
    {
        fRequiredSignalState = Signal::kTime;
    }

    HighPassFilterFFTGenerator::~HighPassFilterFFTGenerator()
    {
    }

    bool HighPassFilterFFTGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        if( aParam->has( "threshold" ) )
        {
        SetThreshold( aParam->get_value< double >( "threshold", fThreshold ) ); 
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
        double CutoffFreq = fThreshold;
        int nwindows = 10;
        int windowsize = aSignal->TimeSize()/nwindows;


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
                    SignalComplex[index][0] = aSignal->SignalTimeComplex()[ ch*aSignal->TimeSize() + nwin*windowsize + index ][0];
                    SignalComplex[index][1] = aSignal->SignalTimeComplex()[ ch*aSignal->TimeSize() + nwin*windowsize + index ][1];
                }

                fftw_execute(ForwardPlan);

                // High Pass FilterFFT
                for( unsigned index = 0; index < windowsize; ++index )
                {
                    if ((index < CutoffFreq/(fAcquisitionRate*1.e6)*windowsize)||(index>0.5*windowsize))
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
                    aSignal->SignalTimeComplex()[ ch*aSignal->TimeSize() + nwin*windowsize + index ][0] = SignalComplex[index][0]/norm;
                    aSignal->SignalTimeComplex()[ ch*aSignal->TimeSize() + nwin*windowsize + index ][1] = SignalComplex[index][1]/norm;
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

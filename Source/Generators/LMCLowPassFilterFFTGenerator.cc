/*
 * LMCLowPassFilterFFTGenerator.cc
 *
 *  Created on: Sept 9, 2016
 *      Author: plslocum after nsoblath
 */

#include "LMCLowPassFilterFFTGenerator.hh"

#include "logger.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "LowPassFilterFFTGenerator" );

    MT_REGISTER_GENERATOR(LowPassFilterFFTGenerator, "lpf-fft");

    LowPassFilterFFTGenerator::LowPassFilterFFTGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &LowPassFilterFFTGenerator::DoGenerateTime ),
        fThreshold( 0.85 ) // fraction of Nyquist above which signals are suppressed.
    {
        fRequiredSignalState = Signal::kTime;
    }

    LowPassFilterFFTGenerator::~LowPassFilterFFTGenerator()
    {
    }

    bool LowPassFilterFFTGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "threshold" ) )
        {
        	SetThreshold( aParam.get_value< double >( "threshold", fThreshold ) );
        }

        return true;
    }

    double LowPassFilterFFTGenerator::GetThreshold() const
    {
        return fThreshold;
    }

    void LowPassFilterFFTGenerator::SetThreshold( double aThreshold )
    {
        fThreshold = aThreshold;
        if (fThreshold > 0.95)
        {
            LERROR(lmclog,"LPF threshold is too close to Nyquist frequency.  It should be < 0.9 .\n");
        	exit(-1);
        }
        return;
    }

// Count leading zeroes in window.
    int LowPassFilterFFTGenerator::GetStartingMargin( Signal* aSignal, int windowsize, int nwin, int ch )
    {
    	int startingMargin = 0;
    	for (int i=0; i<windowsize; i++)
    	{
    		if (fabs(aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + nwin*windowsize + i ][0]) > 0.)
    		{
    			startingMargin = i;
    			break;
    		}
    	}
    	return startingMargin;
    }

// Count trailing zeroes in window.
    int LowPassFilterFFTGenerator::GetEndingMargin( Signal* aSignal, int windowsize, int nwin, int ch )
    {
        int endingMargin = 0;
   	    for (int i=0; i<windowsize; i++)
        {
            if (((nwin+1)*windowsize -i) < 0)
            {
                if (fabs(aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + (nwin+1)*windowsize - i ][0]) > 0.)
                {
                    endingMargin = i;
                    break;
                }
            }
        }
        return endingMargin;
    }



    void LowPassFilterFFTGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    bool LowPassFilterFFTGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool LowPassFilterFFTGenerator::DoGenerateTime( Signal* aSignal )
    {

        LWARN(lmclog,"Calculating LPF-FFT ...\n");

        const unsigned nchannels = fNChannels;

        double CutoffFreq = fThreshold * fAcquisitionRate/2. * 1.e6; // Hz
        double FastNyquist = fAcquisitionRate/2. * 1.e6 * aSignal->DecimationFactor();
        int nwindows = 80;
        int windowsize = aSignal->DecimationFactor()*aSignal->TimeSize()/nwindows;

        for (int ch=0; ch<nchannels; ch++)
        {
            for (int nwin = 0; nwin < nwindows; nwin++)
            {
                int startingMargin = startingMargin = GetStartingMargin(aSignal, windowsize, nwin, ch);
                int endingMargin = GetEndingMargin(aSignal, windowsize, nwin, ch);
                int variableWindowSize = windowsize - startingMargin - endingMargin;

                if (variableWindowSize > 1)
                {
                    fftw_complex *SignalComplex;
                    SignalComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * variableWindowSize );
                    fftw_complex *FFTComplex;
                    FFTComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * variableWindowSize );

                    fftw_plan ForwardPlan;
                    ForwardPlan = fftw_plan_dft_1d(variableWindowSize, SignalComplex, FFTComplex, FFTW_FORWARD, FFTW_ESTIMATE);
                    fftw_plan ReversePlan;
                    ReversePlan = fftw_plan_dft_1d(variableWindowSize, FFTComplex, SignalComplex, FFTW_BACKWARD, FFTW_ESTIMATE);

                    // Construct complex voltage.
                    for( unsigned index = 0; index < variableWindowSize; ++index )
                    {
                        SignalComplex[index][0] = aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + nwin*windowsize + index + startingMargin ][0];
                        SignalComplex[index][1] = aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + nwin*windowsize + index + startingMargin ][1];
                    }

                    fftw_execute(ForwardPlan);

                    // Low Pass FilterFFT
                    for( unsigned index = 0; index < variableWindowSize; ++index )
                    {
                	    if ( (index > variableWindowSize/2.*CutoffFreq/FastNyquist) && (index < variableWindowSize/2. * (1. + (FastNyquist-CutoffFreq)/FastNyquist)))
                        {
                		    FFTComplex[index][0] = 0.;
                		    FFTComplex[index][1] = 0.;
                        }
                    }

                    fftw_execute(ReversePlan);

                    double norm = (double)(variableWindowSize);

                    for( unsigned index = 0; index < variableWindowSize; ++index )
                    {
                        // normalize
                        aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + nwin*windowsize + index + startingMargin ][0] = SignalComplex[index][0]/norm;
                        aSignal->LongSignalTimeComplex()[ ch*aSignal->TimeSize()*aSignal->DecimationFactor() + nwin*windowsize + index + startingMargin ][1] = SignalComplex[index][1]/norm;
                    }

                    fftw_destroy_plan(ForwardPlan);
                    fftw_destroy_plan(ReversePlan);
                    fftw_free( SignalComplex );
                    fftw_free( FFTComplex );
                    SignalComplex = NULL;
                    FFTComplex = NULL;
                }
                else
                {
                    LWARN(lmclog,"LPF-FFT variable window size is too small.\n");
                }

            }  // nwin
        }  // NCHANNELS

        return true;
    }

    bool LowPassFilterFFTGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */

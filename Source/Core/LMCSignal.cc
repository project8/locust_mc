/*
 * LMCSignal.cc
 *
 *  Created on: Feb 3, 2014
 *      Author: nsoblath
 */

#include "LMCSignal.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "Signal" );

    Signal::Signal() :
            fState( kNotInitialized ),
            fTimeSize( 0 ),
            fFreqSize( 0 ),
            fFreqSizeNew( 0 ),
            fDigitalSize( 0 ),
            fSignalTime( NULL ),
            fSignalFreq( NULL ),
            fSignalTimeNew( NULL ),
            fSignalFreqNew( NULL ),
            fLongSignalTimeNew( NULL ),
            fLongSignalFreqNew( NULL ),
            fSignalDigital(),
            fDigitalIsSigned( false ),
            fPlanToFreq(),
            fPlanToTime(),
            fPlanToFreqNew(),
            fPlanToTimeNew()

    {
        fSignalDigital.fUnsigned = NULL;
    }

    Signal::~Signal()
    {
        Reset();
    }

    bool Signal::Initialize( unsigned aTimeSize, unsigned aFFTFlags )
    {
        Reset();

        fTimeSize = aTimeSize;
        fFreqSize = fTimeSize / 2 + 1;
        fFreqSizeNew = fTimeSize;

        fSignalTime = new double [fTimeSize];
        fSignalFreq = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * fFreqSize );
        fSignalTimeNew = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * fTimeSize );
        fSignalFreqNew = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * fFreqSizeNew );
        fLongSignalTimeNew = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * fTimeSize*10 );
        fLongSignalFreqNew = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * fFreqSizeNew*10 );

        ResetValues();

        fPlanToFreq = fftw_plan_dft_r2c_1d( fTimeSize, fSignalTime, fSignalFreq, aFFTFlags);
        fPlanToTime = fftw_plan_dft_c2r_1d( fTimeSize, fSignalFreq, fSignalTime, aFFTFlags);
        fPlanToFreqNew = fftw_plan_dft_1d( fTimeSize, fSignalTimeNew, fSignalFreqNew, FFTW_FORWARD, aFFTFlags);
        fPlanToTimeNew = fftw_plan_dft_1d( fTimeSize, fSignalFreqNew, fSignalTimeNew, FFTW_BACKWARD, aFFTFlags);
        fLongPlanToFreqNew = fftw_plan_dft_1d( fTimeSize*10, fLongSignalTimeNew, fLongSignalFreqNew, FFTW_FORWARD, aFFTFlags);
        fLongPlanToTimeNew = fftw_plan_dft_1d( fTimeSize*10, fLongSignalFreqNew, fLongSignalTimeNew, FFTW_BACKWARD, aFFTFlags);


        fState = kTime;  // pls confused here.
//        fState = kFreq;

        LDEBUG( lmclog, "Signal initialized; time size = " << fTimeSize << "; state = " << fState );

        return true;
    }

    void Signal::Reset()
    {
        delete [] fSignalTime;
        fftw_free( fSignalFreq );
        fftw_free( fSignalTimeNew );
        fftw_free( fSignalFreqNew );
        fftw_free( fLongSignalTimeNew );
        fftw_free( fLongSignalFreqNew );

        if( fDigitalIsSigned ) delete [] fSignalDigital.fSigned;
        else delete [] fSignalDigital.fUnsigned;

        fSignalTime = NULL;
        fSignalFreq = NULL;
        fSignalTimeNew = NULL;
        fSignalFreqNew = NULL;
        fLongSignalTimeNew = NULL;
        fLongSignalFreqNew = NULL;
        fSignalDigital.fUnsigned = NULL;
        fDigitalIsSigned = false;

        fftw_destroy_plan( fPlanToFreq );
        fftw_destroy_plan( fPlanToTime );
        fftw_destroy_plan( fPlanToFreqNew );
        fftw_destroy_plan( fPlanToTimeNew );
        fftw_destroy_plan( fLongPlanToFreqNew );
        fftw_destroy_plan( fLongPlanToTimeNew );


        fTimeSize = 0;
        fFreqSize = 0;
        fFreqSizeNew = 0;
        fDigitalSize = 0;

        fState = kNotInitialized;

        return;
    }

    void Signal::ResetValues()
    {
        if( fSignalTime != NULL )
        {
            for( unsigned index = 0; index < fTimeSize; ++index )
            {
                fSignalTime[index] = 0.;
                fSignalTimeNew[index][0] = 0.;
                fSignalTimeNew[index][1] = 0.;
            }
        }

        if( fSignalFreq != NULL )
        {
            for( unsigned index = 0; index < fFreqSize; ++index )
            {
                fSignalFreq[index][0] = 0.;
                fSignalFreq[index][1] = 0.;
            }
        }


        if( fSignalFreqNew != NULL )
        {
            for( unsigned index = 0; index < fFreqSizeNew; ++index )
            {
                fSignalFreqNew[index][0] = 0.;
                fSignalFreqNew[index][1] = 0.;
            }
        }

        if( fLongSignalTimeNew != NULL )
        {
            for( unsigned index = 0; index < fTimeSize*10; ++index )
            {
                fLongSignalTimeNew[index][0] = 0.;
                fLongSignalTimeNew[index][1] = 0.;
            }
        }

        if( fLongSignalFreqNew != NULL )
        {
            for( unsigned index = 0; index < fFreqSizeNew*10; ++index )
            {
                fLongSignalFreqNew[index][0] = 0.;
                fLongSignalFreqNew[index][1] = 0.;
            }
        }




        if( fDigitalIsSigned )
        {
            if( fSignalDigital.fSigned != NULL )
            {
                for( unsigned index = 0; index < fDigitalSize; ++index )
                {
                    fSignalDigital.fSigned[index] = 0;
                }
            }
        }
        else
        {
            if( fSignalDigital.fUnsigned != NULL )
            {
                for( unsigned index = 0; index < fDigitalSize; ++index )
                {
                    fSignalDigital.fUnsigned[index] = 0;
                }
            }
        }

        return;
    }

    Signal& Signal::Add( const Signal& aSignal, double aScale )
    {
        ToState( fState );

        if( fState == kTime ) AddTime( aSignal, aScale );
        if( fState == kFreq ) AddFreq( aSignal, aScale );

        return *this;
    }
    Signal& Signal::Subtract( const Signal& aSignal, double aScale )
    {
        ToState( fState );

        if( fState == kTime ) SubtractTime( aSignal, aScale );
        if( fState == kFreq ) SubtractFreq( aSignal, aScale );

        return *this;
    }
    Signal& Signal::Multiply( const Signal& aSignal )
    {
        ToState( fState );

        if( fState == kTime ) MultiplyTime( aSignal );
        if( fState == kFreq ) MultiplyFreq( aSignal );

        return *this;
    }
    Signal& Signal::Divide( const Signal& aSignal )
    {
        ToState( fState );

        if( fState == kTime ) DivideTime( aSignal );
        if( fState == kFreq ) DivideFreq( aSignal );

        return *this;
    }

    void Signal::AddTime( const Signal& aSignal, double aScale )
    {
        for( unsigned index = 0; index < fTimeSize; ++index )
        {
            fSignalTime[index] += aSignal.fSignalTime[index] * aScale;
        }
        return;
    }
    void Signal::SubtractTime( const Signal& aSignal, double aScale )
    {
        for( unsigned index = 0; index < fTimeSize; ++index )
        {
            fSignalTime[index] -= aSignal.fSignalTime[index] * aScale;
        }
        return;
    }
    void Signal::MultiplyTime( const Signal& aSignal )
    {
        for( unsigned index = 0; index < fTimeSize; ++index )
        {
            fSignalTime[index] *= aSignal.fSignalTime[index];
        }
        return;
    }
    void Signal::DivideTime( const Signal& aSignal )
    {
        for( unsigned index = 0; index < fTimeSize; ++index )
        {
            fSignalTime[index] /= aSignal.fSignalTime[index];
        }
        return;
    }

    void Signal::AddFreq( const Signal& aSignal, double aScale )
    {
        for( unsigned index = 0; index < fFreqSize; ++index )
        {
            fSignalFreq[index][0] += aSignal.fSignalFreq[index][0] * aScale;
            fSignalFreq[index][1] += aSignal.fSignalFreq[index][1] * aScale;
        }
        return;
    }
    void Signal::SubtractFreq( const Signal& aSignal, double aScale )
    {
        for( unsigned index = 0; index < fFreqSize; ++index )
        {
            fSignalFreq[index][0] -= aSignal.fSignalFreq[index][0] * aScale;
            fSignalFreq[index][1] -= aSignal.fSignalFreq[index][1] * aScale;
        }
        return;
    }
    void Signal::MultiplyFreq( const Signal& aSignal )
    {
        double tempReal, tempImag;
        for( unsigned index = 0; index < fFreqSize; ++index )
        {
            tempReal = fSignalFreq[index][0] * aSignal.fSignalFreq[index][0] - fSignalFreq[index][1] * aSignal.fSignalFreq[index][1];
            tempImag = fSignalFreq[index][0] * aSignal.fSignalFreq[index][1] + fSignalFreq[index][1] * aSignal.fSignalFreq[index][0];
            fSignalFreq[index][0] = tempReal;
            fSignalFreq[index][1] = tempImag;
        }
        return;
    }
    void Signal::DivideFreq( const Signal& aSignal )
    {
        double denom, tempReal, tempImag;
        for( unsigned index = 0; index < fFreqSize; ++index )
        {
            denom = 1. / (aSignal.fSignalFreq[index][0] * aSignal.fSignalFreq[index][0] + aSignal.fSignalFreq[index][1] * aSignal.fSignalFreq[index][1]);
            tempReal = fSignalFreq[index][0] * aSignal.fSignalFreq[index][0] + fSignalFreq[index][1] * aSignal.fSignalFreq[index][1];
            tempImag = fSignalFreq[index][1] * aSignal.fSignalFreq[index][0] - fSignalFreq[index][0] * aSignal.fSignalFreq[index][1];
            fSignalFreq[index][0] = tempReal * denom;
            fSignalFreq[index][1] = tempImag * denom;
        }
        return;
    }

    unsigned Signal::TimeSize() const
    {
        return fTimeSize;
    }

    unsigned Signal::FreqSize() const
    {
        return fFreqSize;
    }

    unsigned Signal::DigitalSize() const
    {
        return fDigitalSize;
    }

    bool Signal::ToState( Signal::State aState )
    {
        if( fState == kNotInitialized )
        {
            LERROR( lmclog, "Signal is not initialized" );
            return false;
        }
        if( fState == aState ) return true;
        if( aState == kTime ) return FFTToTime();
        if( aState == kFreq ) return FFTToFreq();
        if( aState == kDigital )
        {
            LERROR( lmclog, "Please use Signal::ToDigital to convert to a digital signal" );
            return false;
        }
        LERROR( lmclog, "Unknown state requested: " << aState );
        return false;
    }

    Signal::State Signal::GetState() const
    {
        return fState;
    }

    bool Signal::ToTime()
    {
        return ToState( kTime );
    }

    bool Signal::ToFreq()
    {
        return ToState( kFreq );
    }

//    bool Signal::ToDigital( uint64_t* anArray, unsigned aDigSize )    
    bool Signal::ToDigital( int8_t* anArray, unsigned aDigSize )  // pls
    {
        if( fDigitalIsSigned ) delete [] fSignalDigital.fSigned;
        else delete [] fSignalDigital.fUnsigned;
        fSignalDigital.fSigned = anArray;
        fDigitalSize = aDigSize;
        fState = kDigital;
        fDigitalIsSigned = true;
        return true;
    }

    bool Signal::ToDigital( uint8_t* anArray, unsigned aDigSize )  // pls
    {
        if( fDigitalIsSigned ) delete [] fSignalDigital.fSigned;
        else delete [] fSignalDigital.fUnsigned;
        fSignalDigital.fUnsigned = anArray;
        fDigitalSize = aDigSize;
        fState = kDigital;
        fDigitalIsSigned = false;
        return true;
    }

    bool Signal::FFTToTime()
    {
        LDEBUG( lmclog, "Performing reverse FFT to the time domain" );
        fftw_execute( fPlanToTime );
        fState = kTime;
        return true;
    }

    bool Signal::FFTToFreq()
    {
        LDEBUG( lmclog, "Performing forward FFT to frequency domain" );
        fftw_execute( fPlanToFreq );
        fState = kFreq;
        return true;  // pls
    }

    const double* Signal::SignalTime() const
    {
        return fSignalTime;
    }

    double* Signal::SignalTime()
    {
        return fSignalTime;
    }

    const fftw_complex* Signal::SignalTimeNew() const
    {
        return fSignalTimeNew;
    }

    fftw_complex* Signal::SignalTimeNew()
    {
        return fSignalTimeNew;
    }



    const fftw_complex* Signal::LongSignalTimeNew() const
    {
        return fLongSignalTimeNew;
    }

    fftw_complex* Signal::LongSignalTimeNew()
    {
        return fLongSignalTimeNew;
    }


/*
    double Signal::SignalTime( unsigned anIndex ) const
    {
        return fSignalTime[ anIndex ];
    }

    double& Signal::SignalTime( unsigned anIndex )
    {
        return fSignalTime[ anIndex ];
    }
*/
    const fftw_complex* Signal::SignalFreq() const
    {
        return fSignalFreq;
    }

    fftw_complex* Signal::SignalFreq()
    {
        return fSignalFreq;
    }
/*
    const fftw_complex& Signal::SignalFreq( unsigned anIndex ) const
    {
        return fSignalFreq[ anIndex ];
    }

    fftw_complex& Signal::SignalFreq( unsigned anIndex )
    {
        return fSignalFreq[ anIndex ];
    }
*/
//    const uint64_t* Signal::SignalDigital() const    
    const int8_t* Signal::SignalDigitalS() const  // pls
    {
        return fSignalDigital.fSigned;
    }

    const uint8_t* Signal::SignalDigitalUS() const  // pls
    {
        return fSignalDigital.fUnsigned;
    }

//    uint64_t* Signal::SignalDigital()    
    int8_t* Signal::SignalDigitalS()  // pls
    {
        return fSignalDigital.fSigned;
    }

    uint8_t* Signal::SignalDigitalUS()  // pls
    {
        return fSignalDigital.fUnsigned;
    }
/*
//    uint64_t Signal::SignalDigital( unsigned anIndex ) const    
    uint8_t Signal::SignalDigital( unsigned anIndex ) const  // pls
    {
        return fSignalDigital[ anIndex ];
    }

//    uint64_t& Signal::SignalDigital( unsigned anIndex )    
    uint8_t& Signal::SignalDigital( unsigned anIndex )  // pls
    {
        return fSignalDigital[ anIndex ];
    }
*/
} /* namespace locust */

/*
 * LMCSignal.cc
 *
 *  Created on: Feb 3, 2014
 *      Author: nsoblath
 */

#include "LMCSignal.hh"

#include "LMCLogger.hh"

namespace locust
{
    LMCLOGGER( lmclog, "Signal" );

    Signal::Signal() :
            fState( kNotInitialized ),
            fTimeSize( 0 ),
            fFreqSize( 0 ),
            fSignalTime( NULL ),
            fSignalFreq( NULL ),
            fPlanToFreq(),
            fPlanToTime()
    {
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

        fSignalTime = new double [fTimeSize];
        fSignalFreq = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * fFreqSize );

        ResetValues();

        fPlanToFreq = fftw_plan_dft_r2c_1d( fTimeSize, fSignalTime, fSignalFreq, aFFTFlags);
        fPlanToTime = fftw_plan_dft_c2r_1d( fTimeSize, fSignalFreq, fSignalTime, aFFTFlags);

        fState = kTime;

        return true;
    }

    void Signal::Reset()
    {
        delete [] fSignalTime;
        fftw_free( fSignalFreq );

        fftw_destroy_plan( fPlanToFreq );
        fftw_destroy_plan( fPlanToTime );

        fTimeSize = 0;
        fFreqSize = 0;

        fState = kNotInitialized;

        return;
    }

    void Signal::ResetValues()
    {
        for( unsigned index = 0; index < fTimeSize; ++index )
        {
            fSignalTime[index] = 0.;
        }

        for( unsigned index = 0; index < fFreqSize; ++index )
        {
            fSignalFreq[index][0] = 0.;
            fSignalFreq[index][1] = 0.;
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

    bool Signal::ToState( Signal::State aState )
    {
        if( fState == kNotInitialized )
        {
            LMCERROR( lmclog, "Signal is not initialized" );
            return false;
        }
        if( fState == aState ) return true;
        if( aState == kTime ) return ToTime();
        if( aState == kFreq ) return ToFreq();
        LMCERROR( lmclog, "Unknown state requested: " << aState );
        return false;
    }

    Signal::State Signal::GetState() const
    {
        return fState;
    }

    bool Signal::ToTime()
    {
        fftw_execute( fPlanToTime );
        return true;
    }

    bool Signal::ToFreq()
    {
        fftw_execute( fPlanToFreq );
        return false;
    }

    const double* Signal::SignalTime() const
    {
        return fSignalTime;
    }

    double* Signal::SignalTime()
    {
        return fSignalTime;
    }

    double Signal::SignalTime( unsigned anIndex ) const
    {
        return fSignalTime[ anIndex ];
    }

    double& Signal::SignalTime( unsigned anIndex )
    {
        return fSignalTime[ anIndex ];
    }

    const fftw_complex* Signal::SignalFreq() const
    {
        return fSignalFreq;
    }

    fftw_complex* Signal::SignalFreq()
    {
        return fSignalFreq;
    }

    const fftw_complex& Signal::SignalFreq( unsigned anIndex ) const
    {
        return fSignalFreq[ anIndex ];
    }

    fftw_complex& Signal::SignalFreq( unsigned anIndex )
    {
        return fSignalFreq[ anIndex ];
    }

} /* namespace locust */

/*
 * LMCSignal.hh
 *
 *  Created on: Feb 3, 2014
 *      Author: nsoblath
 */

#ifndef LMCSIGNAL_HH_
#define LMCSIGNAL_HH_

#include <complex.h>
#include <fftw3.h>

namespace locust
{
    class DigitalSignalCore;

    class Signal
    {
        public:
            enum State
            {
                kNotInitialized,
                kTime,
                kFreq,
                kDigital
            };

        public:
            Signal();
            virtual ~Signal();

            bool Initialize( unsigned aTimeSize, unsigned aFFTFlags = FFTW_ESTIMATE );
            void Reset();
            void ResetValues();

            Signal& Add( const Signal& aSignal, double aScale = 1. );
            Signal& Subtract( const Signal& aSignal, double aScale = 1. );
            Signal& Multiply( const Signal& aSignal );
            Signal& Divide( const Signal& aSignal );

            unsigned TimeSize() const;
            unsigned FreqSize() const;
            unsigned DigitalSize() const;

            bool ToState( State aState );
            State GetState() const;

            bool ToTime();
            bool ToFreq();
            bool ToDigital( uint64_t* anArray, unsigned aDigSize );

            const double* SignalTime() const;
            double* SignalTime();

            double SignalTime( unsigned anIndex ) const;
            double& SignalTime( unsigned anIndex );

            const fftw_complex* SignalFreq() const;
            fftw_complex* SignalFreq();

            const fftw_complex& SignalFreq( unsigned anIndex ) const;
            fftw_complex& SignalFreq( unsigned anIndex );

            const uint64_t* SignalDigital() const;
            uint64_t* SignalDigital();

            uint64_t SignalDigital( unsigned anIndex ) const;
            uint64_t& SignalDigital( unsigned anIndex );

        private:
            bool FFTToTime();
            bool FFTToFreq();

            void AddTime( const Signal& aSignal, double aScale );
            void SubtractTime( const Signal& aSignal, double aScale );
            void MultiplyTime( const Signal& aSignal );
            void DivideTime( const Signal& aSignal );

            void AddFreq( const Signal& aSignal, double aScale );
            void SubtractFreq( const Signal& aSignal, double aScale );
            void MultiplyFreq( const Signal& aSignal );
            void DivideFreq( const Signal& aSignal );

            State fState;

            unsigned fTimeSize;
            unsigned fFreqSize;
            unsigned fDigitalSize;

            double* fSignalTime;
            fftw_complex* fSignalFreq;
            uint64_t* fSignalDigital;

            fftw_plan fPlanToFreq;
            fftw_plan fPlanToTime;

    };

} /* namespace locust */

#endif /* LMCSIGNAL_HH_ */

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

    class Signal
    {
        public:
            enum State
            {
                kNotInitialized,
                kTime,
                kFreq
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

            bool ToState( State aState );
            State GetState() const;

            const double* SignalTime() const;
            double* SignalTime();

            double SignalTime( unsigned anIndex ) const;
            double& SignalTime( unsigned anIndex );

            const fftw_complex* SignalFreq() const;
            fftw_complex* SignalFreq();

            const fftw_complex& SignalFreq( unsigned anIndex ) const;
            fftw_complex& SignalFreq( unsigned anIndex );

        private:
            bool ToTime();
            bool ToFreq();

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

            double* fSignalTime;
            fftw_complex* fSignalFreq;

            fftw_plan fPlanToFreq;
            fftw_plan fPlanToTime;

    };

} /* namespace locust */

#endif /* LMCSIGNAL_HH_ */

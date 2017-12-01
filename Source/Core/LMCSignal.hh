/*
 * LMCSignal.hh
 *
 *  Created on: Feb 3, 2014
 *      Author: nsoblath
 */

#ifndef LMCSIGNAL_HH_
#define LMCSIGNAL_HH_

#include <fftw3.h>
#include <complex.h>
#include <stdint.h>


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
//            bool ToDigital( uint64_t* anArray, unsigned aDigSize );            
            bool ToDigital( int8_t* anArray, unsigned aDigSize );
            bool ToDigital( uint8_t* anArray, unsigned aDigSize );  // pls

            const double* SignalTime() const;
            double* SignalTime();

            //double SignalTime( unsigned anIndex ) const;
            //double& SignalTime( unsigned anIndex );

            const fftw_complex* SignalTimeComplex() const;
            fftw_complex* SignalTimeComplex();

            const fftw_complex* LongSignalTimeComplex() const;
            fftw_complex* LongSignalTimeComplex();

            const fftw_complex* LongSignalFreqNew() const;
            fftw_complex* LongSignalFreqNew();

            const fftw_complex* SignalFreq() const;
            fftw_complex* SignalFreq();

            const fftw_complex* SignalFreqNew() const;
            fftw_complex* SignalFreqNew();


            //const fftw_complex& SignalFreq( unsigned anIndex ) const;
            //fftw_complex& SignalFreq( unsigned anIndex );

//            const uint64_t* SignalDigital() const;
            const int8_t* SignalDigitalS() const;
            const uint8_t* SignalDigitalUS() const;  // pls
//            uint64_t* SignalDigital();
            int8_t* SignalDigitalS();
            uint8_t* SignalDigitalUS();  // pls

            bool GetDigitalIsSigned() const;

//            uint64_t SignalDigital( unsigned anIndex ) const;            
            //uint8_t SignalDigital( unsigned anIndex ) const;  // pls
//            uint64_t& SignalDigital( unsigned anIndex );            
            //uint8_t& SignalDigital( unsigned anIndex );  // pls

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
            unsigned fFreqSizeNew;
            unsigned fDigitalSize;

            double* fSignalTime;
            fftw_complex* fSignalFreq;
            fftw_complex* fSignalTimeComplex;
            fftw_complex* fSignalFreqNew;
            fftw_complex* fLongSignalTimeComplex;
            fftw_complex* fLongSignalFreqNew;


            //            uint64_t* fSignalDigital;
            union SignalDigital
            {
                int8_t* fSigned;
                uint8_t* fUnsigned;
            } fSignalDigital;
            bool fDigitalIsSigned;

            fftw_plan fPlanToFreq;
            fftw_plan fPlanToTime;
            fftw_plan fPlanToFreqNew;
            fftw_plan fPlanToTimeComplex;
            fftw_plan fLongPlanToFreqNew;
            fftw_plan fLongPlanToTimeComplex;

    };

    inline bool Signal::GetDigitalIsSigned() const
    {
        return fDigitalIsSigned;
    }

} /* namespace locust */

#endif /* LMCSIGNAL_HH_ */

/*
 * LMCRunLengthCalculator.hh
 *
 *  Created on: Feb 7, 2014
 *      Author: nsoblath
 */

#ifndef LMCRUNLENGTHCALCULATOR_HH_
#define LMCRUNLENGTHCALCULATOR_HH_

#include "LMCVisitor.hh"

namespace scarab
{
    class param_node;
}

namespace locust
{

    /*!
     @class RunLengthCalculator
     @author N. S. Oblath

     @brief Calculates the length of the run according to direct user requests and generator settings

     @details

     Available configuration options:
     - "n-records": unsigned -- Set the run length by the number of records
     - "duration": double -- Set the run length by total time in seconds
     - "record-size": unsigned -- Number of elements in a record
     - "acquisition-rate": double -- Set the acquisition rate in MHz (also sets bin width)
     - "bin-width" double -- Set the bin width rate in ns (also sets acquisition rate)
    */
    class RunLengthCalculator : public GeneratorVisitor
    {
        public:
            enum RunLengthState
            {
                kDecided = -1,
                kByRecords = 1,
                kByDuration = 2,
                kByGenerators = 3,
                kUnknown = 100
            };

        public:
            RunLengthCalculator();
            virtual ~RunLengthCalculator();

            bool Configure( const scarab::param_node& aNode );

            bool VisitGenerators();

            bool CalculateRunLength();

            RunLengthState GetState() const;
            /// Set the run length state iff state < fState
            void SetRunLengthState( RunLengthState state );
            /// Set the run length state regardless of the current state
            void OverrideRunLengthState( RunLengthState state );

            unsigned GetNRecords() const;
            void SetNRecords( unsigned recs );

            double GetDuration() const;
            void SetDuration( double duration );

            unsigned GetRecordSize() const;
            void SetRecordSize( unsigned size );

            unsigned GetSampleSize() const;
            void SetSampleSize( unsigned size );

            unsigned GetNChannels() const;
            void SetNChannels( unsigned size );

            void SetFirstGenerator( const Generator* firstGen );

            RunLengthState GetByGeneratorsState() const;
            unsigned GetByGeneratorsNRecords() const;
            double GetByGeneratorsDuration() const;

            double GetAcquisitionRate() const;
            /// sets the acquisition rate and the bin width (= 1/ar)
            void SetAcquisitionRate( double ar );

            double GetBinWidth() const;
            /// sets the bin width and the acquisition rate (= 1/bw)
            void SetBinWidth( double bw );

        private:
            void Visit( const KassSignalGenerator* );
            void Visit( const FreeFieldSignalGenerator* );
            void Visit( const PatchSignalGenerator* );
            void Visit( const GaussianNoiseGenerator* );
            void Visit( const FakeTrackSignalGenerator* );
            void Visit( const PlaneWaveSignalGenerator* );
            void Visit( const AntennaSignalGenerator* );
            void Visit( const TestSignalGenerator* );
            void Visit( const LowPassFilterFFTGenerator* );
            void Visit( const HighPassFilterFFTGenerator* );
            void Visit( const LocalOscillatorGenerator* );
            void Visit( const DecimateSignalGenerator* );
            void Visit( const DipoleSignalGenerator* );
            void Visit( const TurnstileSignalGenerator* );
            void Visit( const Digitizer* );

            RunLengthState fState;

            unsigned fNRecords;
            unsigned frlcChannels;


            double fDuration; // sec

            const Generator* fFirstGenerator;
            RunLengthState fByGeneratorsState;
            unsigned fByGeneratorsNRecords;
            double fByGeneratorsDuration; // sec

            unsigned fRecordSize;
            unsigned fSampleSize;
            double fBinWidth; // sec
            double frlcAcquisitionRate; // MHz

    };

} /* namespace locust */

#endif /* LMCRUNLENGTHCALCULATOR_HH_ */

/*
 * LMCRunLengthCalculator.hh
 *
 *  Created on: Feb 7, 2014
 *      Author: nsoblath
 */

#ifndef LMCRUNLENGTHCALCULATOR_HH_
#define LMCRUNLENGTHCALCULATOR_HH_

#include "LMCVisitor.hh"

namespace locust
{
    class ParamNode;

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

            bool Configure( const ParamNode* aNode );

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
            void Visit( const GaussianNoiseGenerator* );

            RunLengthState fState;

            unsigned fNRecords;

            double fDuration; // sec

            const Generator* fFirstGenerator;
            RunLengthState fByGeneratorsState;
            unsigned fByGeneratorsNRecords;
            double fByGeneratorsDuration; // sec

            unsigned fRecordSize;
            double fBinWidth; // sec
            double fAcquisitionRate; // MHz

    };

} /* namespace locust */

#endif /* LMCRUNLENGTHCALCULATOR_HH_ */

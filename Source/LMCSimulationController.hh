/*
 * LMCSimulationController.hh
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#ifndef LMCSIMULATIONCONTROLLER_HH_
#define LMCSIMULATIONCONTROLLER_HH_

namespace locust
{
    class Generator;
    class ParamNode;

    class SimulationController
    {
        public:
            SimulationController();
            virtual ~SimulationController();

            void Configure( const ParamNode* aNode );

            unsigned GetNRecords() const;
            /// Set the number of records; also modifies the duration
            void SetNRecords( unsigned recs );

            double GetDuration() const;
            /// Set the duration: rounds down to an integer number of records, and sets the duration accordingly
            void SetDuration( double duration );

            unsigned GetRecordSize() const;
            /// Set the record size; also modifies the duration
            void SetRecordSize( unsigned size );

            double GetBinWidth() const;
            /// Set the bin width; also modifies the duration
            void SetBinWidth( double bw );

            void SetFirstGenerator( const Generator* firstGen );

            void Run() const;

        private:
            unsigned fNRecords;
            double fDuration; // sec

            unsigned fRecordSize;
            double fBinWidth; // sec

            const Generator* fFirstGenerator;
    };

} /* namespace locust */

#endif /* LMCSIMULATIONCONTROLLER_HH_ */

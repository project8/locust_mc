/*
 * LMCEggWriter.hh
 *
 *  Created on: Feb 6, 2014
 *      Author: nsoblath
 */

#ifndef LMCEGGWRITER_HH_
#define LMCEGGWRITER_HH_

#include "M3Monarch.hh"

#include "member_variables.hh"

namespace scarab
{
    class param_node;
}

namespace locust
{
    class Signal;

    class RunLengthCalculator;
    class Digitizer;

    /*!
     @class EggWriter
     @author N. S. Oblath

     @brief Uses the Monarch library to write simulated signals to an egg file

     @details

     Available configuration options:
     - "filename": string -- filename of the output egg file
     - "date": string -- date/time string
     - "description": string -- informative run description
    */
    class EggWriter
    {
        public:
            enum State
            {
                kClosed,
                kPrepared,
                kWriting,
            };

        public:
            EggWriter();
            virtual ~EggWriter();

            bool Configure( const scarab::param_node& aNode );

            mv_referrable( std::string, filename );
            mv_referrable( std::string, date );
            mv_referrable( std::string, description );

        public:
            bool PrepareEgg( const RunLengthCalculator* aRLC, const Digitizer* aDig );

            bool WriteRecord( const Signal* aSignal, bool is_new_acq );

            bool FinalizeEgg();

            void IncrementAcquisitionId();

            mv_accessible_noset( monarch3::TimeType, record_length );
            mv_accessible_noset( monarch3::RecordIdType, record_id );
            mv_accessible_noset( monarch3::TimeType, record_time );
            mv_accessible_noset( uint64_t, record_n_bytes );

        private:
            monarch3::Monarch3* f_monarch;
            monarch3::M3Stream* f_stream;
            monarch3::M3Record* f_record;

        private:
            State f_state;

    };

} /* namespace locust */

#endif /* LMCEGGWRITER_HH_ */

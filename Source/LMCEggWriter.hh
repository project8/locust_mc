/*
 * LMCEggWriter.hh
 *
 *  Created on: Feb 6, 2014
 *      Author: nsoblath
 */

#ifndef LMCEGGWRITER_HH_
#define LMCEGGWRITER_HH_

#include "Monarch.hpp"

namespace locust
{
    class ParamNode;
    class Signal;

    /*!
     @class EggWriter
     @author N. S. Oblath

     @brief Uses the Monarch library to write simulated signals to an egg file

     @details

     Available configuration options:
     - "filename": string -- filename of the output egg file
     - "date": string -- date/time string
     - "description": string -- informative run description
     - "run-type": unsigned -- set the run type
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

            bool Configure( const ParamNode* aNode );

            const std::string& GetFilename() const;
            void SetFilename( const std::string& filename );

            const std::string& GetDate() const;
            void SetDate( const std::string& date );

            const std::string& GetDescription() const;
            void SetDescription( const std::string& desc );

            monarch::RunType GetRunType() const;
            void SetRunType( monarch::RunType runType );

            unsigned GetBitDepth() const;
            bool SetBitDepth( unsigned bitDepth );

            unsigned GetDataTypeSize() const;

            double GetVoltageMin() const;
            void SetVoltageMin( double vMin );

            double GetVoltageRange() const;
            void SetVoltageRange( double vRange );

            double GetAcquisitionRate() const;
            void SetAcquisitionRate( double rate );

            double GetDuration() const;
            void SetDuration( double duration );

            unsigned GetRecordSize() const;
            void SetRecordSize( unsigned size );

        private:
            State fState;

            // This info will be configured directly
            std::string fFilename;
            std::string fDate;
            std::string fDescription;
            monarch::RunType fRunType;

            // This info will come from the digitizer
            unsigned fBitDepth;
            unsigned fDataTypeSize;
            double fVoltageMin;
            double fVoltageRange;

            // This info will come from the run-length calculator
            double fAcquisitionRate; // MHz
            double fDuration; // s
            unsigned fRecordSize;

            // Calculated at runtime
            monarch::TimeType fRecordLength; // ns

        public:
            bool PrepareEgg();

            bool WriteRecord( const Signal* aSignal );

            bool FinalizeEgg();

            monarch::Monarch* GetMonarch() const;

            void IncrementAcquisitionId();

        private:
            monarch::Monarch* fMonarch;

            monarch::MonarchRecordBytes* fRecord;

            monarch::AcquisitionIdType fAcquisitionId;
            monarch::RecordIdType fRecordCounter;
            monarch::TimeType fRecordTime;

            uint64_t fRecordNBytes;

    };

} /* namespace locust */

#endif /* LMCEGGWRITER_HH_ */
